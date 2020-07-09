You can estimate the raw number of variant in each sample using the following command:

```bash
for i in SVvariants/*bcf; do
  basename $i .bcf
  bcftools view --no-header $i | wc -l
done | paste - - | column -t
```

You should get:

|NA12878|123|
|NA12891|155|
|NA12892|132|

