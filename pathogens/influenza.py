import csv
import dataclasses
import datetime

from pathogen_properties import *
from populations import us_population

background = """Influenza is a respiratory virus.  While we call it one thing,
it has four variants (A/B/C/D) that don't form a clade.  Only A, B, and C
infect humans, and C is rare enough it's not usually counted."""

FLU_A = TaxID(11320)
FLU_B = TaxID(11520)

pathogen_chars = PathogenChars(
    na_type=NAType.RNA,
    enveloped=Enveloped.ENVELOPED,
    taxids=frozenset((FLU_A, FLU_B)),
    subtaxids=frozenset((FLU_A, FLU_B)),
    names_by_taxid={
        FLU_A: "Influenza A",
        FLU_B: "Influenza B",
    },
    selection=SelectionRound.ROUND_1,
)


# https://ndc.services.cdc.gov/wp-content/uploads/2021/02/MMWR_Week_overview.pdf
# "The first day of any MMWR week is Sunday. MMWR week numbering is
#  sequential beginning with 1 and incrementing with each week to a
#  maximum of 52 or 53. MMWR week #1 of an MMWR year is the first week
#  of the year that has at least four days in the calendar year. For
#  example, if January 1 occurs on a Sunday, Monday, Tuesday or
#  Wednesday, the calendar week that includes January 1 would be MMWR
#  week #1. If January 1 occurs on a Thursday, Friday, or Saturday, the
#  calendar week that includes January 1 would be the last MMWR week of
#  the previous year (#52 or #53). Because of this rule, December 29,
#  30, and 31 could potentially fall into MMWR week #1 of the following
#  MMWR year."
# Another way of saying this is that in years where Jan 1 falls on Thursday,
# Friday, Saturday, or Sunday then the first MMWR week of the year starts on
# the first Sunday of the year, while if it falls on Monday, Tuesday, or
# Wednesday then the first week starts with the last Sunday of the previous
# year.
def parse_mmwr_week(year: int, week: int) -> datetime.date:
    first_day_of_year = datetime.date(year, 1, 1)
    first_sunday_of_year = datetime.date(
        year, 1, 7 - first_day_of_year.weekday()
    )
    assert first_sunday_of_year.weekday() == 6

    if first_day_of_year.weekday() in [3, 4, 5, 6]:
        # First sunday of the year is included in MMWR week 1.
        return first_sunday_of_year + datetime.timedelta(weeks=week - 1)
    else:
        # First sunday of the year is not included in MMWR week 1.
        return first_sunday_of_year + datetime.timedelta(weeks=week - 2)


WeeklyData = dict[str, dict[datetime.date, tuple[int, int]]]


def load_weekly_data() -> WeeklyData:
    # State -> Week -> (PositiveA, PositiveB)
    output: WeeklyData = {}

    # Downloaded 2024-08-15 from
    # https://gis.cdc.gov/grasp/fluview/fluportaldashboard.html with options:
    #   Select Data Source:
    #      [x] WHO/NREVSS   [ ] ILINet
    #   [x] State
    #   Select Regions: [x] Select All
    #   Select Seasons:
    #     [x] 2023-24
    #     [x] 2022-23
    #     [x] 2021-22
    #     [x] 2020-21
    #     [x] 2019-20
    #
    # This is CDC data on the number of positive tests by week, for Flu A vs
    # Flu B.
    #
    # This comes from clinical labs; I can't find week-level data for the
    # public health labs, but there's already a lot here.
    with open(
        prevalence_data_filename("CDC_WHO_NREVSS_Clinical_Labs.csv")
    ) as inf:
        cols = None
        for row in csv.reader(inf):
            if row[0].startswith("*"):
                continue  # initial comment

            if not cols:
                cols = row
                continue

            region = row[cols.index("REGION")]
            year = int(row[cols.index("YEAR")])
            mmwr_week = int(row[cols.index("WEEK")])
            parsed_start = parse_mmwr_week(year, mmwr_week)
            total_tests = row[cols.index("TOTAL SPECIMENS")]
            positive_a = row[cols.index("TOTAL A")]
            positive_b = row[cols.index("TOTAL B")]

            if total_tests == positive_a == positive_b == "X":
                continue

            if region not in output:
                output[region] = {}

            output[region][parsed_start] = (
                int(positive_a),
                int(positive_b),
            )
    return output


infections_2019_2020 = IncidenceAbsolute(
    annual_infections=36_000_000,
    confidence_interval=(28_663_846, 71_197_342),
    coverage_probability=0.95,  # 95% UI
    start_date="2019-07-01",
    end_date="2020-07-01",
    tag="us-2019-2020",
    country="United States",
    # "The overall burden of influenza (flu) for the 2019-2020 was an estimated
    # 36 million flu-related illnesses"
    source="https://www.cdc.gov/flu/about/burden/2019-2020.html#:~:text=The%20overall%20burden%20of%20influenza%20(flu)%20for%20the%202019%2D2020%20was%20an%20estimated%C2%A036%20million%C2%A0flu%2Drelated%20illnesses",
)

infections_2021_2022 = IncidenceAbsolute(
    annual_infections=9_000_000,
    confidence_interval=(7_683_440, 15_541_059),
    coverage_probability=0.95,  # 95% UI
    start_date="2021-07-01",
    end_date="2022-07-01",
    tag="us-2021-2022",
    country="United States",
    # "The overall burden of influenza (flu) for the 2021-2022 season was an
    # estimated 9 million flu illnesses"
    source="https://www.cdc.gov/flu/about/burden/2021-2022.htm#:~:text=The%20overall%20burden%20of%20influenza%20(flu)%20for%20the%202021%2D2022%20season%20was%20an%20estimated%C2%A09%20million%C2%A0flu%20illnesses",
)


def compare_incidence_to_positive_tests(
    official_incidence: IncidenceAbsolute, weekly_data: WeeklyData
) -> Scalar:
    annual_infections = 0
    min_date = None
    max_date = None
    for state in weekly_data:
        for parsed_start in weekly_data[state]:
            parsed_end = parsed_start + datetime.timedelta(weeks=1)
            assert official_incidence.parsed_start
            assert official_incidence.parsed_end
            if (
                parsed_start >= official_incidence.parsed_start
                and parsed_end <= official_incidence.parsed_end
            ):
                if not min_date or min_date > parsed_start:
                    min_date = parsed_start

                if not max_date or max_date < parsed_end:
                    max_date = parsed_end

                positive_a, positive_b = weekly_data[state][parsed_start]
                annual_infections += positive_a + positive_b

    assert min_date
    assert max_date
    clinical_labs_infections = IncidenceAbsolute(
        annual_infections=annual_infections,
        tag="us-%s-%s" % (min_date.year, max_date.year),
    )

    return official_incidence / clinical_labs_infections


def estimate_incidences() -> list[IncidenceRate]:
    # State -> Week -> (PositiveA, PositiveB)
    weekly_data = load_weekly_data()

    # We can't just go from positive tests to incidence because many people
    # will get flu but never be tested, or their test won't make it back to the
    # CDC.  For each flu season we use the overall CDC incidence estimate and
    # our count of positive tests to get an estimate of underreporting.

    underreporting_2019_2020 = compare_incidence_to_positive_tests(
        infections_2019_2020, weekly_data
    )
    underreporting_2021_2022 = compare_incidence_to_positive_tests(
        infections_2021_2022, weekly_data
    )

    # The CDC didn't estimate an annual infections for 2020-2021, but assume
    # underreporting is the average of the two years we do have.
    underreporting_2020_2021 = Scalar.average(
        underreporting_2019_2020, underreporting_2021_2022
    )

    def get_underreporting(
        start: datetime.date, end: datetime.date
    ) -> Optional[Scalar]:
        assert infections_2019_2020.parsed_start
        assert infections_2019_2020.parsed_end
        assert infections_2021_2022.parsed_start
        assert infections_2021_2022.parsed_end

        if (
            start >= infections_2019_2020.parsed_start
            and start <= infections_2019_2020.parsed_end
        ):
            return underreporting_2019_2020
        elif (
            start >= infections_2021_2022.parsed_start
            and start <= infections_2021_2022.parsed_end
        ):
            return underreporting_2021_2022
        elif (
            start >= infections_2019_2020.parsed_end
            and start <= infections_2021_2022.parsed_start
        ):
            return underreporting_2020_2021
        else:
            return None

    incidences = []
    for state in ["California", "Ohio"]:
        for parsed_start in weekly_data[state]:
            positive_a, positive_b = weekly_data[state][parsed_start]
            parsed_end = parsed_start + datetime.timedelta(weeks=1)

            for taxid, weekly_count in [
                (FLU_A, positive_a),
                (FLU_B, positive_b),
            ]:
                adjusted_weekly_count = float(weekly_count)
                if parsed_start.year <= 2019:
                    continue

                underreporting = get_underreporting(parsed_start, parsed_end)
                if not underreporting:
                    continue

                is_pseudocount = False
                if weekly_count == 0:
                    # See comment on QUANTITY_WHEN_NONE_OBSERVED.
                    adjusted_weekly_count = QUANTITY_WHEN_NONE_OBSERVED
                    is_pseudocount = True

                incidence = IncidenceAbsolute(
                    annual_infections=adjusted_weekly_count * 52,
                    country="United States",
                    state=state,
                    date=parsed_start.isoformat(),
                )

                incidences.append(
                    dataclasses.replace(
                        (
                            incidence.to_rate(
                                us_population(
                                    state=state, year=parsed_start.year
                                )
                            )
                            * underreporting
                        ),
                        taxid=taxid,
                        is_pseudocount=is_pseudocount,
                    )
                )

    return incidences


def estimate_prevalences() -> list[Prevalence]:
    return []
