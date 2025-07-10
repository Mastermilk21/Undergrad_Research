from astroquery.gaia import Gaia
import pandas as pd
import time

df = pd.read_csv("matched_transient_coords_no_duplicates.csv")  # Replace with actual path

results = []

for index, row in df.iterrows():
    name = row["name"]
    ra = row["RA"]
    dec = row["Dec"]
    
    query = f"""
    SELECT TOP 1
        source_id,
        ra, dec,
        parallax, parallax_error,
        phot_g_mean_mag,
        phot_bp_mean_mag,
        phot_rp_mean_mag,
        DISTANCE(POINT('ICRS', ra, dec), POINT('ICRS', {ra}, {dec})) * 3600 AS angular_distance_arcsec
    FROM gaiadr3.gaia_source
    WHERE 1=CONTAINS(
        POINT('ICRS', ra, dec),
        CIRCLE('ICRS', {ra}, {dec}, {2/3600.0})
    )
    ORDER BY angular_distance_arcsec ASC
    """

    try:
        job = Gaia.launch_job_async(query)
        r = job.get_results()

        if len(r) > 0:
            best = r[0]
            results.append({
                "name": name,
                "RA": ra,
                "Dec": dec,
                "source_id": best["source_id"],
                "parallax": best["parallax"],
                "parallax_error": best["parallax_error"],
                "phot_g_mean_mag": best["phot_g_mean_mag"],
                "phot_bp_mean_mag": best["phot_bp_mean_mag"],
                "phot_rp_mean_mag": best["phot_rp_mean_mag"],
                "angular_distance_arcsec": best["angular_distance_arcsec"],
            })
        else:
            results.append({
                "name": name,
                "RA": ra,
                "Dec": dec,
                "source_id": None,
                "parallax": None,
                "parallax_error": None,
                "phot_g_mean_mag": None,
                "phot_bp_mean_mag": None,
                "phot_rp_mean_mag": None,
                "angular_distance_arcsec": None,
            })
    except Exception as e:
        print(f"Error on {name}: {e}")
        continue

    time.sleep(0.5)  # be gentle on the server

# Write to file
output_df = pd.DataFrame(results)
output_df.to_csv("gaia_matches.csv", index=False)
print("✅ Gaia query complete. Results saved to 'gaia_matches.csv'.")
