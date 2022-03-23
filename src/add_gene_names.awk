NR==FNR{
        if($3=="gene"){
                gene_names[$10] = $18
        }
        next
}

$3!="gene"{
        print($0, "gene_name " gene_names[$10])
        next
}
