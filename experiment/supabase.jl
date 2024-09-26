using HTTP
using JSON

struct SupabaseClient
    url::String
    key::String
end

function get_table(sc::SupabaseClient, tablename::String)
    response = HTTP.get(
        "$(sc.url)/rest/v1/$tablename",
        headers = [
            "apikey" => sc.key,
            "Authorization" => "Bearer $(sc.key)",
            "Content-Type" => "application/json"
        ]
    )

    return JSON.parse(String(response.body))
end

function update_heights(client, heights, tablename)
    table = get_table(client, tablename)
    dict = table[1]
    for i in 1:10
        dict[string(i)] += heights[i]
    end
    dict["inf"] += heights[11]

    updatedPayload = JSON.json(dict)

    patchUrl = "$(client.url)/rest/v1/$tablename?id=eq.1"
    patchResponse = HTTP.patch(
        patchUrl,
        headers = [
            "apikey" => client.key,
            "Authorization" => "Bearer $(client.key)",
            "Content-Type" => "application/json",
            "Prefer" => "return=representation"
        ],
        body = updatedPayload
    )

    # println("PATCH Status: ", patchResponse.status)
    return nothing
end

function add_polys(client, all_dicts, polystablename)
    postUrl = "$(client.url)/rest/v1/$polystablename"
    # println("all_dicts: $all_dicts")
    payload = JSON.json(all_dicts)
    # println("payload: $payload")
    response = HTTP.post(
        postUrl,
        headers = [
            "apikey" => client.key,
            "Authorization" => "Bearer $(client.key)",
            "Content-Type" => "application/json",
            "Prefer" => "return=representation"
        ],
        body = payload
    )

    # println("POST Status: ", response.status)
end
