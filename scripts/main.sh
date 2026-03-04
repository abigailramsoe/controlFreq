name="SMDB_$(date +"%Y%m%d_%H%M%S")"
rm tmpdir/*
curl -L -o "tmpdir/${name}.zip" -X POST -d 'checkbox_smdb=on' http://dandyweb01fl.unicph.domain:5100/download_merged_standardized
unzip -o "tmpdir/${name}.zip" -d tmpdir/
tsv=$(find tmpdir -type f -name "*.tsv" | head -n 1)
mv "$tsv" "smdb/$name.tsv"
rm tmpdir/*
echo "SMDB file downloaded and controls extracted: smdb/$name.tsv"

echo "Getting control data and generating report..."
Rscript scripts/getControls.R $(realpath "smdb/$name.tsv")

