import pandas as pd
from astroquery.heasarc import Heasarc
from astropy.coordinates import SkyCoord
import astropy.units as u

heasarc = Heasarc()

xray_catalogs = {
    'ROSAT All-Sky Survey': 'rosmaster',
    'Chandra Source Catalog': 'csc',
    'XMM-Newton SSC': 'xmmssc',
    # The correct catalog name for Swift XRT (you may try 'swiftmastr' or check HEASARC docs)
    'Swift XRT Point Source Catalog': 'swiftmastr'  
}

input_csv = 'matched_transient_coords_no_duplicates.csv'  # Replace with your actual file
output_csv = 'xray_crossmatch_results.csv'

df = pd.read_csv(input_csv)

results = []

for idx, row in df.iterrows():
    name = row['name']
    ra = row['RA']
    dec = row['Dec']
    coord = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    print(f"Checking {name} at RA={ra}, Dec={dec}")

    for cat_name, cat_code in xray_catalogs.items():
        try:
            query_result = heasarc.query_region(coord, catalog=cat_code, radius=2*u.arcsec)
            if query_result is not None and len(query_result) > 0:
                for source in query_result:
                    results.append({
                        'transient_name': name,
                        'catalog': cat_name,
                        'source_name': source.get('Name', ''),
                        'ra_source': source['RA'],
                        'dec_source': source['DEC'],
                        'flux': source.get('FLUX', ''),
                        'flux_error': source.get('FLUX_ERR', ''),
                        'other_info': str(source)
                    })
        except Exception as e:
            print(f"Error querying {cat_name} for {name}: {e}")

if results:
    df_results = pd.DataFrame(results)
    df_results.to_csv(output_csv, index=False)
    print(f"Cross-match complete. Results saved to {output_csv}")
else:
    print("No X-ray matches found for any sources.")
