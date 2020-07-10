You can estimate the raw number of variant in each sample using the following command:

```bash
for i in SVvariants/*bcf; do
  basename $i .bcf
  bcftools view --no-header $i | wc -l
done | paste - - | column -t
```

You should get:

|Sample|count|
|--|--|
|NA12878|212|
|NA12891|227|
|NA12892|207|

