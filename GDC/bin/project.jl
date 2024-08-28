using GDC
using HTTP
using Debugger

let endpoint = "projects",
    filter = (Field("project.project_id") == "TCGA-LUAD"),
    fields = ["file_id"]
    global body = GDC.GDCString( filter, fields )
    println.( [filter, fields, body] );
    response = HTTP.post( "https://api.gdc.cancer.gov/$endpoint", ["Content-Type" => "application/json"], body )
end

body = "{\n\t\"filters\":{\n\t\t\"op\":\"=\",\n\t\t\"content\":{\n\t\t\t\"field\":\"project.project_id\",\n\t\t\t\"value\":\"TCGA-LUAD\"\n\t\t}\n\t},\n\t\"fields\":\"file_id\"\n}"

HTTP.request( "POST", "https://api.gdc.cancer.gov/projects", ["Content-Type" => "application/json"], body )
sleep(2)
HTTP.request( "POST", "https://api.gdc.cancer.gov/projects", ["Content-Type" => "application/json"], body )
