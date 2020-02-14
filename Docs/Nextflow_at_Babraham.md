# Nextflow Pipelines at the Babraham Institute

#### General comments: 


## Pipelines:

## Single Progam Pipelines:


All preconfigured pipelines take one additional argument, which has to be exactly in the following form to work:

```
--toolname_args "'--additional_option value --extra_flag etc.'"
```

So as an example, you could run specific trimming in Trim Galore like so:

```
--trim_galore_args "'--clip_r1 10 --clip_r2 10 --nextera'"
```
