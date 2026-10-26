nextflow.enable.dsl=2

// Groovy port of the row-grouping logic in bin/houseley_group_setup_multirun_structure.py.
// This MUST stay in lock-step with that script's read_run_info()/format_cell()/check_if_dual()
// functions, since the directory names computed here have to match the directory names the
// Python script actually creates on disk.

def parseHouseleySamplesheet(String csvPath, boolean noHeader = false) {

    def groups = [:]  // dirName -> [runType: ..., genome: ..., groupId: ...]

    def lines = new File(csvPath).readLines()
    def dataLines = noHeader ? lines : lines.drop(1)

    dataLines.each { rawLine ->
        if (rawLine.trim().isEmpty()) {
            return
        }

        def cells = parseHouseleyCsvLine(rawLine)
        // mirror: line = [format_cell(cell) for cell in row if format_cell(cell)]
        def line = cells.collect { houseleyFormatCell(it) }.findAll { it }

        if (line.size() < 3) {
            return
        }

        def dual = houseleyIsDualBarcode(line[2])

        def groupId
        def runType
        def genome
        def dirName

        // dirName is bounded to groupID/run type/genome only -- columns after genome (e.g.
        // TrAEL indexes, names) are for other purposes and must not end up in the directory name.
        // Must stay in sync with the same bound in houseley_group_setup_multirun_structure.py.
        if (dual) {
            groupId = line[3]
            runType = line.size() > 4 ? line[4] : null
            genome  = line.size() > 5 ? line[5] : null
            dirName = line[3..Math.min(5, line.size() - 1)].join('_')
        } else {
            groupId = line[2]
            runType = line.size() > 3 ? line[3] : null
            genome  = line.size() > 4 ? line[4] : null
            dirName = line[2..Math.min(4, line.size() - 1)].join('_')
        }

        // Columns 7/8 (raw, 1-indexed) are the TrAEL index list (e.g. "5,6,7") and the sample
        // names for those indexes (e.g. "Alg70bp1, 70 in ESAG8, TbCas9"), in the same order.
        // Read from the raw (unformatted) cells, not the reformatted 'line' array: houseleyFormatCell()
        // would hyphen-join the index list's commas (e.g. "5-6-7") rather than leaving the plain
        // concatenated digits nf_traelseq's --index expects (e.g. "567"), and would similarly mangle
        // the names list, which needs its commas intact for houseley_rename_trael_indexes.py. These
        // columns' positions in the raw CSV are fixed regardless of dual/non-dual barcoding, since
        // the sample sheet always has the same number of columns per row (an unused Barcode 2
        // column is still present, just empty, for non-dual rows).
        def rawIndexCell = cells.size() > 6 ? cells[6] : null
        def indexDigits = rawIndexCell ? rawIndexCell.replaceAll(/[,\s]/, "") : null
        indexDigits = indexDigits ? indexDigits : null  // treat an empty result as absent
        def rawNamesCell = cells.size() > 7 ? cells[7] : null

        if (!groups.containsKey(dirName)) {
            groups[dirName] = [runType: runType, genome: genome, groupId: groupId,
                                traelIndexDigits: [] as Set, traelRenames: []]
        }

        // Multiple sample sheet rows can share the same dirName (e.g. several TrAEL samples
        // pooled into one group folder), each with its own barcode/index/names -- so this is
        // accumulated for every row, not just the first one seen for a dirName like runType/genome
        // above (which are the same for every row sharing a dirName, by construction of dirName
        // itself).
        if (indexDigits) {
            indexDigits.each { d -> groups[dirName].traelIndexDigits.add(d as String) }
            groups[dirName].traelRenames.add([barcode: line[1], index: indexDigits, names: rawNamesCell])
        }
    }

    // --index for nf_traelseq needs every index digit used by any sample sharing the group's
    // folder (nf_traelseq just keeps files matching any of them; per-sample renaming afterwards
    // is scoped separately by barcode via traelRenames).
    groups.each { dirName, info ->
        info.traelIndex = info.traelIndexDigits ? info.traelIndexDigits.sort().join("") : null
    }

    return groups
}

// RFC4180-ish single CSV line parser (handles quoted fields containing commas, and "" as an
// escaped quote within a quoted field) -- matches Python's csv.reader closely enough for the
// Houseley sample sheets, which only use quoting to protect commas inside the sample name.
def parseHouseleyCsvLine(String line) {
    def fields = []
    def sb = new StringBuilder()
    boolean inQuotes = false
    boolean skipNext = false
    def chars = line.toCharArray()

    chars.eachWithIndex { c, idx ->
        if (skipNext) {
            skipNext = false
        } else if (inQuotes) {
            if (c == ('"' as char)) {
                if (idx + 1 < chars.length && chars[idx + 1] == ('"' as char)) {
                    sb.append('"')
                    skipNext = true
                } else {
                    inQuotes = false
                }
            } else {
                sb.append(c)
            }
        } else {
            if (c == ('"' as char)) {
                inQuotes = true
            } else if (c == (',' as char)) {
                fields.add(sb.toString())
                sb = new StringBuilder()
            } else {
                sb.append(c)
            }
        }
    }
    fields.add(sb.toString())
    return fields
}

def houseleyFormatCell(String cell) {
    cell = cell.trim()
    if (cell.contains(',')) {
        def items = cell.split(',').collect { it.trim() }
        cell = items.join('-')
    }
    cell = cell.replace(' ', '_')
    cell = cell.replace('(', '').replace(')', '')
    return cell
}

def houseleyIsDualBarcode(String s) {
    return s ==~ /^[CATG]*$/
}

// True if the group's fastq folder contains an R2 file (i.e. the raw data is paired-end).
// Trim Galore preserves R1/R2 pairing into _val_1/_val_2, so checking the raw files here is
// valid whether or not the group goes through the Vazyme pre-trim step first.
def houseleyIsPairedEnd(File dirPath) {
    if (! dirPath.exists()) { return false }
    def fastqFiles = dirPath.list().findAll { it.endsWith(".fastq.gz") }
    return fastqFiles.any { it =~ /_R2\.fastq\.gz$/ }
}

// Names of the group folder's fastq.gz files that contain no reads at all. A sample with zero
// reads after demultiplexing still gets a (20 byte) fastq.gz file, so file size alone isn't a
// reliable test -- decompress just far enough to see whether there's a single byte of content.
// A file that can't be read at all isn't counted as empty: let the dispatched pipeline fail on
// it loudly instead.
def houseleyEmptyFastqs(File dirPath) {
    if (! dirPath.exists()) { return [] }
    def fastqFiles = dirPath.listFiles().findAll { it.name.endsWith(".fastq.gz") }
    return fastqFiles.findAll { f ->
        try {
            new java.util.zip.GZIPInputStream(new FileInputStream(f)).withCloseable { it.read() == -1 }
        } catch (IOException e) {
            false
        }
    }.collect { it.name }.sort()
}
