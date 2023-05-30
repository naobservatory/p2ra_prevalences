from pathogen_properties import *
from populations import us_population

background = """Hepatitis B is a liver infection caused by the Hepatitis B
 virus. It is transmitted through birth and contact with infected blood or
 bodily fluids. Hepatitis B infection increases the risk for hepatocellular
 carcinoma"""

# TODO:
# 1. Add in-house NHANES estimate. Specifically, we could present NHANES data
# for HBV-core antibody (evidence of past infection), and HBV-surface antigen
# (evidence of active acute or chronic infection):
# https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/HEPBD_J.htm
# 2. Incorporate an estimate based on mortality date. This data can be acquired
# from the CDC WONDER database (https://wonder.cdc.gov/mcd-icd10-provisional.html).
# Parameters for relevant query is:
# - Select time period of death": 2022, and
# - Select underlying cause of death: B16 (Acute hepatitis B).
# All other fields can be left as default.
# 3. Consider adding an underreporting factor for chronic infections. If I don't
#  do this we might drop reported chronic infections, as these do not correspond
# to prevalence by a factor of around 100x.


pathogen_chars = PathogenChars(
    na_type=NAType.DNA,
    enveloped=Enveloped.ENVELOPED,
    taxid=TaxID(10407),
)

us_population_2019 = Population(
    people=328.2 * 1e6,
    country="United States",
    date="2019",
    source="https://data.census.gov/table?q=us+population+2019&tid=ACSDP1Y2019.DP05",
)

cdc_acute_underreporting_factor_2020 = Scalar(
    scalar=6.5,
    confidence_interval=(3.7, 15.9),
    coverage_probability=0.95,
    country="United States",
    source="https://www.cdc.gov/hepatitis/statistics/2020surveillance/introduction/technical-notes.htm#:~:text=each%20reported%20case%20of%20acute%20hepatitis%20B%20represents%206.5%20estimated%20infections%20(95%25%20bootstrap%20CI%3A%203.7%E2%80%9315.9",
)


cdc_estimated_acute_2019 = IncidenceAbsolute(
    annual_infections=20_700,
    # During 2019, a total of 3,192 acute hepatitis B cases were reported to
    # CDC, resulting in 20,700 estimated infections (95% CI: 11,800–50,800)
    # after adjusting for case underascertainment and underreporting.
    confidence_interval=(11_800, 50_800),  # 95% Bootstrap Confidence Interval
    coverage_probability=0.95,
    country="United States",
    # I picked 2019 rates, as estimated rates in 2020 were 30% lower, even
    # while deaths stayed the same, pointing to undercounting due to
    # COVID-19: https://www.cdc.gov/hepatitis/statistics/2020surveillance/introduction/national-profile.htm#:~:text=hepatitis%20B%20transmission.-,Data%20from%20death%20certificates%20filed%20in%20the%20vital%20records%20offices%20of,same%20as%20the%20rate%20during%202019%20(0.42%20deaths%20per%20100%2C000%20population).,-Hepatitis%20C
    date="2019",
    source="https://www.cdc.gov/hepatitis/statistics/2019surveillance/Introduction.htm#Technical:~:text=20%2C700%20estimated%20infections%20(95%25%20CI%3A%2011%2C800%E2%80%9350%2C800)",
    methods="https://www.cdc.gov/hepatitis/statistics/2019surveillance/Introduction.htm#Technical:~:text=To%20account%20for,CI%3A%2011.0%E2%80%9347.4).",
)

cdc_reported_chronic_2019 = IncidenceRate(
    annual_infections_per_100k=5.9,  # no CI provided
    # For chronic incidence we do not currently use an underreporting factor.
    country="United States",
    date="2019",
    source="https://www.cdc.gov/hepatitis/statistics/2019surveillance/Introduction.htm#ref09:~:text=A%20total%20of%2013%2C859%20new%20cases%20of%20chronic%20hepatitis%20B%20were%20reported%20to%20CDC%20during%202019%2C%20corresponding%20to%20a%20rate%20of%205.9%20cases%20per%20100%2C000%20population.",
)


estimated_chronic_prevalence_us_2020 = PrevalenceAbsolute(
    infections=1_721_027,
    # This is their mid estimate in table 4, second to last row. The estimate
    # listed in their paper and on the last line of the table doesn't add up.
    # They say they added 303'237 (point estimate US-born) + 1'076'069 (lower
    # estimate Foreign-born) + 72'013 (point estimate Non-NHANES). But that
    # adds up to 1'450'320, which isn't the number they give (1.59M). We are
    # contacting the authors on this.
    confidence_interval=(1_249_055, 2_491_435),  # These numbers represent a
    # "low" and "high" estimate, labels that aren't further defined. Numbers
    # can be found in table 4, second to last row.
    active=Active.LATENT,
    country="United States",
    date="2020",
    source="https://journals.lww.com/ajg/fulltext/2020/09000/prevalence_of_chronic_hepatitis_b_virus_infection.20.aspx#:~:text=Table%204.%3A%20Estimated%20prevalence%20of%20chronic%20HBV%20in%20the%20United%20Statesa",
    methods="https://journals.lww.com/ajg/fulltext/2020/09000/prevalence_of_chronic_hepatitis_b_virus_infection.20.aspx#:~:text=Panel%20members%20researched,in%20the%20US.",
)


ohio_reported_acute_incidence_2021 = IncidenceRate(
    # This source reports both acute and total (acute+chronic) reported incidence.
    # We report both acute and chronic incidence seperately.
    annual_infections_per_100k=0.9,
    country="United States",
    state="Ohio",
    date="2021",
    # Note that this source has yearly data for 2017-2021 if it's ever needed
    source="https://odh.ohio.gov/know-our-programs/viral-hepatitis/data-statistics/hbv_5yr_report",
)

ohio_reported_chronic_incidence_2021 = IncidenceRate(
    # Same source as above, here we subtract the acute incidence from
    # the total incidence.
    # For chronic incidence we do not currently use an underreporting factor.
    annual_infections_per_100k=16.5 - 0.9,
    # Importantly, these are newly identified cases and do not represent prevalence.
    country="United States",
    county="Ohio",
    date="2021",
    source="https://odh.ohio.gov/know-our-programs/viral-hepatitis/data-statistics/hbv_5yr_report",
)


def estimate_prevalences() -> list[Prevalence]:
    return [
        estimated_chronic_prevalence_us_2020.to_rate(us_population(year=2020))
    ]


def estimate_incidences() -> list[IncidenceRate]:
    return [
        cdc_estimated_acute_2019.to_rate(us_population_2019),
        cdc_reported_chronic_2019,
        ohio_reported_acute_incidence_2021
        * cdc_acute_underreporting_factor_2020,
        ohio_reported_chronic_incidence_2021,
        # We do not have an underreporting factor for chronic incidence
    ]
