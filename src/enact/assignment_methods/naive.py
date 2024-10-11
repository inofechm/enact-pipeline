# Naive method. Only using the bins unique to each cell (overlapping bins omitted)


def naive_assignment(result_spatial_join):
    # Naive method. Only using the bins unique to each cell (overlapping bins omitted)
    result_spatial_join = result_spatial_join[result_spatial_join["unique_bin"]]
    result_spatial_join["weight"] = 1
    return result_spatial_join
