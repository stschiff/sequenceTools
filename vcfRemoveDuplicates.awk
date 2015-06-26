$0 ~ /^#/ {print}
$0 !~ /^#/ {
    if(!($1 == last_chrom && $2 == last_pos))
        print
    last_chrom = $1
    last_pos = $2
}