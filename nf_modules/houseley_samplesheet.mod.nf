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

        if (dual) {
            dirName = line[3..-1].join('_')
            groupId = line[3]
            runType = line.size() > 4 ? line[4] : null
            genome  = line.size() > 5 ? line[5] : null
        } else {
            dirName = line[2..-1].join('_')
            groupId = line[2]
            runType = line.size() > 3 ? line[3] : null
            genome  = line.size() > 4 ? line[4] : null
        }

        if (!groups.containsKey(dirName)) {
            groups[dirName] = [runType: runType, genome: genome, groupId: groupId]
        }
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
