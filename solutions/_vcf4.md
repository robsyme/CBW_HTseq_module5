You can estimate the raw number of variant by type in each sample using the following command:

```bash
for bcf in SVvariants/*.bcf; do
  basename bcf .bcf
  bcftools view --no-header bcf | awk '{print $5}' | sort | uniq -c
done | paste - - - -
```

|Sample|Deletion|Duplication|Inversion|
|--|--|--|--|
|NA12878|92|15|16|
|NA12891|134|12|9|
|NA12892|113|11|8|
