using HTTP
using CSV
using DataFrames
using GDC

filters = (Field("files.access") == "open") &
    (Field("files.data_category") == "simple nucleotide variation") &
    (Field("files.experimental_strategy") == "WXS")

fields = [
    "file_id",
    "file_name",
    "data_type",
    "data_format",
]

body = GDCString( filters, fields, size=10_000, format="tsv" )

response = HTTP.request( "POST", "https://api.gdc.cancer.gov/files", ["Content-Type" => "application/json"], body );
result = CSV.read( response.body, DataFrame );
