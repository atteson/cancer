using HTTP
using CSV
using DataFrames

body = """
{
    "filters":{
        "op":"and",
        "content":[
            {
                "op":"in",
                "content":{
                    "field":"cases.submitter_id",
                    "value":[
                        "TCGA-CK-4948",
                        "TCGA-D1-A17N",
                        "TCGA-4V-A9QX",
                        "TCGA-4V-A9QM"
                    ]
                }
            },
            {
                "op":"=",
                "content":{
                    "field":"files.data_type",
                    "value":"Gene Expression Quantification"
                }
            }
        ]
    },
    "format":"csv",
    "fields":"file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,analysis.workflow_type,cases.project.project_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id",
    "size":"1000"
}
"""

response = HTTP.request( "POST", "https://api.gdc.cancer.gov/files", ["Content-Type" => "application/json"], body );
result = CSV.read( response.body, DataFrame )
