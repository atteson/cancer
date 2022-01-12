using GDC

filters = (Field("files.access") == "open") &
    (Field("files.data_category") == "simple nucleotide variation") &
    (Field("files.experimental_strategy") == "WXS")

result = find_files( filters, size=10_000, format="tsv" );
