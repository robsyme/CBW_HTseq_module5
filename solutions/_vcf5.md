You can estimate the raw number of variant by type after merging each sample using the following command:

```
for bcf in SVvariants/*.bcf; do
  basename $bcf .bcf
  bcftools view $bcf \
  | awk '/^[^#]/ {print $5}' \
  | sort \
  | uniq -c
done | paste - - - -
```

You should get the following numbers:

|Sample|Deletion|Duplication|Inversion|
|--|--|--|--|
|NA12878|92|15|16|
|NA12891|134|12|9|
|NA12892|113|11|8|
|Merged|86|11|14|
