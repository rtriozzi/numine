import pandas

def import_cafana_data(
    PATH,
    COLUMNS,
    FLUX_MODE,
    **kwargs,
):
    # import data
    df = pandas.read_csv(
        PATH,
        names = COLUMNS,
        **kwargs
    )
    df['nominal_flux_mode'] = FLUX_MODE

    return df