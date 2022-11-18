from pysam import VariantFile

vfc_handler: VariantFile = VariantFile("testing.vcf", mode="rb")
print(vfc_handler)

for variant in vfc_handler:
    print(6 >> 1)
    break
