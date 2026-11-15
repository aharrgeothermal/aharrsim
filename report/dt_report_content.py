"""Static content for the development-types report: what each type is (from the source table), the example that was run,
how it was modelled, and the schematic specification drawn by the page's SVG renderer.

Schematic spec (all depths in km, x in 0..1 across the drawing):
  depth      : depth of the drawing's floor
  ticks      : depth ticks to label (km)
  isotherms  : [(z, label)]           rock temperature labels drawn at the right edge
  bands      : [dict(z, h, label, kind)]   kind: aquifer | store | hot | steam
  wells      : [dict(x, z, kind, dir, label, pump, lateral, n_lat, lat_z)]
               kind: single | coaxial | utube | coaxial_rev ; dir: up | down | loop
  fractures  : dict(x1, x2, z, n, h)   vertical fractures between two wells (h in km)
  units      : [dict(x, label, kind)]  surface plant boxes; kind: orc | flash | steam | hp | hx | heat | sep | heater
  flows      : [dict(a, b, label, hot)] arrows at surface between wellheads ('w0', 'w1', ...) and units ('u0', ...)
  notes      : [dict(x, z, text)]      free labels in the ground
"""

GROUPS = [
    ("Conventional", [1, 2, 3, 3.5]),
    ("Surface waste heat", [4]),
    ("Closed loop, shallow", [5, 6]),
    ("Closed loop, deep", [7, 8, 9]),
    ("EGS / HDR", [10, 11, 12, 13, 14]),
    ("Oil and gas assets", [15, 16, 17]),
    ("Supercritical and ultradeep", [18, 19, 20]),
]

TYPES = {
    1: dict(
        short="Hot sedimentary aquifer doublet",
        sheet="Low-temperature (below about 90 °C) water from a permeable sedimentary aquifer, produced and reinjected by a doublet; "
              "heat for district networks, greenhouses and aquaculture. The source table quotes under 2 MWe and under 5 MWth per well.",
        example="Paris Basin Dogger type: 1.8 km, 70 °C aquifer 20 m thick (porosity 0.15), 68 kg/s produced by a line-shaft pump set at "
                "250 m, returned at 40 °C to a district-heating network, injector 1.2 km away.",
        model="Production and injection wells in <code>mode: production</code> / <code>injection</code> with the heat loss up the 9⅝ in casing; "
              "the aquifer as a Gringarten–Sauty doublet (thermal breakthrough with caprock and bedrock losses), the injector's bottom-hole "
              "temperature fed to the aquifer model and the aquifer temperature fed to the producer at every time step.",
        look="Thermal breakthrough: the produced water starts to cool once the injected front reaches the producer.",
        plant="heat only; binary checked",
        spec=dict(depth=2.2, ticks=[0, 1, 2], isotherms=[(1.8, "70 °C")],
                  bands=[dict(z=1.75, h=0.12, label="Dogger aquifer, 20 m", kind="aquifer")],
                  wells=[dict(x=0.3, z=1.8, kind="single", dir="up", label="producer", pump=0.25),
                         dict(x=0.72, z=1.8, kind="single", dir="down", label="injector")],
                  units=[dict(x=0.51, label="heat network", kind="heat")],
                  flows=[dict(a="w0", b="u0", label="63–69 °C", hot=True), dict(a="u0", b="w1", label="40 °C", hot=False)],
                  notes=[dict(x=0.51, z=1.0, text="1.2 km apart")])),
    2: dict(
        short="Pumped medium-temperature well, binary plant",
        sheet="Medium-temperature (90–150 °C and above) liquid reservoirs produced with downhole pumps to binary (ORC) plants; the Great Basin "
              "type. The source table quotes 2–6 MWe (under 7) and under 15 MWth per well.",
        example="Nevada type: 2.0 km well into a 165 °C convective reservoir declining 0.5 K a year, 75 kg/s lifted by a line-shaft pump at "
                "450 m to a 1.2 MPa wellhead, binary plant returning the brine at 70 °C.",
        model="Single production well with a constant-temperature convective reservoir (<code>ConstantReservoir</code> with decline), "
              "productivity index 40 kg/s/MPa against 18 MPa, downhole pump at 450 m; GEOPHIRES practical binary correlation at 15 °C "
              "ambient, ORC cycle model as a check.",
        look="The wellhead temperature sits a few degrees below the reservoir early on and the gap closes as the casing warms; the pump "
             "takes 10–20 % of the output.",
        plant="binary (GEOPHIRES practical)",
        spec=dict(depth=2.3, ticks=[0, 1, 2], isotherms=[(2.0, "165 °C")],
                  bands=[dict(z=1.9, h=0.3, label="convective reservoir", kind="hot")],
                  wells=[dict(x=0.35, z=2.0, kind="single", dir="up", label="producer", pump=0.45),
                         dict(x=0.75, z=2.0, kind="single", dir="down", label="injector")],
                  units=[dict(x=0.55, label="binary plant", kind="orc")],
                  flows=[dict(a="w0", b="u0", label="147–161 °C", hot=True), dict(a="u0", b="w1", label="70 °C", hot=False)])),
    3: dict(
        short="Self-flowing high-temperature well, flash plant",
        sheet="High-temperature liquid-dominated reservoirs (above 200 °C) in volcanic settings: wells self-flow a two-phase mixture to a "
              "separator and a single- or double-flash plant. The source table quotes 3–25 MWe per well (up to 30–50), 12.5 MWth typical.",
        example="Wairakei type: 1.5 km well with a 250 °C liquid feed (9 MPa reservoir, productivity index 60 kg/s/MPa), self-flowing "
                "40 kg/s to a 0.9 MPa wellhead, flashing in the riser; separator at wellhead pressure.",
        model="Single production well with the drift-flux two-phase riser, feed enthalpy at reservoir pressure, wellhead-pressure boundary; "
              "single and double flash with the Baumann rule, 40 °C condenser, 1 % NCG.",
        look="The wellhead temperature is pinned at saturation for 0.9 MPa (175 °C); the steam fraction (0.16) carries the energy.",
        plant="double flash",
        spec=dict(depth=1.8, ticks=[0, 1], isotherms=[(1.5, "250 °C")],
                  bands=[dict(z=1.35, h=0.4, label="liquid-dominated reservoir", kind="hot")],
                  wells=[dict(x=0.38, z=1.5, kind="single", dir="up", label="self-flowing", flash=0.5),
                         dict(x=0.8, z=1.5, kind="single", dir="down", label="reinjection")],
                  units=[dict(x=0.6, label="flash plant", kind="flash")],
                  flows=[dict(a="w0", b="u0", label="175 °C, x 0.16", hot=True), dict(a="u0", b="w1", label="brine", hot=False)])),
    3.5: dict(
        short="Vapour-dominated field, dry steam",
        sheet="Vapour-dominated reservoirs (The Geysers, Larderello) produce dry or slightly superheated steam straight to a condensing "
              "turbine. The source table's 3–25 MWe per well applies.",
        example="Geysers type: 2.5 km well into a 240 °C saturated-steam reservoir at 3.3 MPa, 15 kg/s of steam to a 0.8 MPa wellhead.",
        model="Single production well fed with saturated-steam enthalpy (productivity index 15 kg/s/MPa), steam-phase riser hydraulics and "
              "casing heat loss; condensing turbine on the wellhead stream (separator removes the 1 % of condensate), 40 °C condenser.",
        look="Heat loss up the casing condenses about 1 % of the steam; the turbine output per kg of steam is five times that of a flash well.",
        plant="steam turbine",
        spec=dict(depth=2.9, ticks=[0, 1, 2], isotherms=[(2.5, "240 °C")],
                  bands=[dict(z=2.3, h=0.5, label="steam reservoir, 3.3 MPa", kind="steam")],
                  wells=[dict(x=0.4, z=2.5, kind="single", dir="up", label="steam well", steam=True),
                         dict(x=0.8, z=2.5, kind="single", dir="down", label="condensate injection")],
                  units=[dict(x=0.6, label="steam turbine", kind="steam")],
                  flows=[dict(a="w0", b="u0", label="170 °C steam", hot=True), dict(a="u0", b="w1", label="condensate", hot=False)])),
    4: dict(
        short="Modular ORC on surface waste heat",
        sheet="Small modular binary units (Climeon type, about 150 kWe each) on low-temperature streams that already exist at surface: "
              "separated brine, engine cooling water, industrial waste heat. The source table quotes 355 kWe per module, 1–2 MWe stacked.",
        example="A 110 °C stream of 50 kg/s at 1 MPa (for instance separated brine), cooled to 70 °C; no wells.",
        model="Plant model only: GEOPHIRES practical binary correlation and the ORC cycle model (R245fa / isobutane) at 15 °C ambient.",
        look="A fixed-temperature source: the output depends only on the stream and the ambient temperature.",
        plant="binary (GEOPHIRES practical)",
        spec=dict(depth=0.0, ticks=[], surface_only=True,
                  units=[dict(x=0.22, label="brine / waste heat", kind="sep"), dict(x=0.58, label="ORC modules", kind="orc"),
                         dict(x=0.86, label="return", kind="heat")],
                  flows=[dict(a="u0", b="u1", label="110 °C, 50 kg/s", hot=True), dict(a="u1", b="u2", label="70 °C", hot=False)])),
    5: dict(
        short="Ground-source heat pump, U-tube borehole",
        sheet="Shallow closed loops (50–300 m) for heat pumps: the borehole supplies a few kilowatts at near-ground temperature. The source "
              "table quotes under 10 kW per unit.",
        example="150 m single U-tube (32 mm SDR11 pipe, 60 mm shank spacing, grout 1.5 W/m/K) in 12 °C ground, 0.3 kg/s at 2 °C inlet for a "
                "182-day heating season, ten seasons; heat pump lifting to 45 °C.",
        model="<code>mode: utube</code>: both legs coupled to the borehole wall through the multipole delta circuit, seasonal schedule, "
              "unit-cell far field at 120 m; heat pump at 45 % of Carnot between the mean loop temperature and 45 °C.",
        look="Season by season the ground cools a little; the loop settles near 5 °C outlet and 29 W/m.",
        plant="heat pump (consumes electricity)",
        spec=dict(depth=0.18, ticks=[0, 0.05, 0.1, 0.15], isotherms=[(0.1, "12 °C")],
                  wells=[dict(x=0.45, z=0.15, kind="utube", dir="loop", label="U-tube")],
                  units=[dict(x=0.72, label="heat pump → 45 °C", kind="hp")],
                  flows=[dict(a="w0", b="u0", label="5 °C", hot=True), dict(a="u0", b="w0", label="2 °C", hot=False, back=True)])),
    6: dict(
        short="Shallow coaxial borehole exchanger",
        sheet="Deep borehole heat exchangers of a few hundred metres (coaxial) feeding heat pumps for buildings. The source table quotes "
              "under 10 kW per installation up to a few MWth for fields of them.",
        example="400 m coaxial exchanger: 6⅝ in casing with a 90 mm PE inner pipe, 2 kg/s at 5 °C inlet, 10 °C surface and 30 °C/km, ten years; "
                "heat pump to 45 °C.",
        model="The coaxial solver with an uninsulated inner pipe (k 0.4 W/m/K), unit-cell far field at 150 m; heat pump as for type 5.",
        look="The uninsulated pipe short-circuits part of the heat; 27 W/m after ten years.",
        plant="heat pump (consumes electricity)",
        spec=dict(depth=0.46, ticks=[0, 0.1, 0.2, 0.3, 0.4], isotherms=[(0.4, "22 °C")],
                  wells=[dict(x=0.45, z=0.4, kind="coaxial", dir="loop", label="coaxial BHE")],
                  units=[dict(x=0.74, label="heat pump → 45 °C", kind="hp")],
                  flows=[dict(a="w0", b="u0", label="6 °C", hot=True), dict(a="u0", b="w0", label="5 °C", hot=False, back=True)])),
    7: dict(
        short="Single-well deep coaxial retrofit in a hot field",
        sheet="A closed coaxial loop run in an existing hot well (GreenFire GreenLoop at Coso): water circulated down the annulus and up an "
              "insulated tubing, no fluid exchange with the reservoir. The source table carries claims of 2–9 MWe per well and a 2 MWe retrofit.",
        example="3.0 km well in an 80 °C/km field (255 °C at TD), 9⅝ in casing to 1.5 km and 7 in below, 3½ × 4½ in vacuum-insulated "
                "tubing (k 0.03 W/m/K), 5 and 10 kg/s injected at 40 °C against 2 MPa.",
        model="The validated coaxial solver (annulus down, tubing up), conduction only, 30 years.",
        look="Even in 255 °C rock a single conduction-fed well gives 0.7 MWth at 10 kg/s and water below 60 °C.",
        plant="binary (GEOPHIRES practical)",
        spec=dict(depth=3.4, ticks=[0, 1, 2, 3], isotherms=[(3.0, "255 °C"), (1.5, "135 °C")],
                  bands=[dict(z=2.6, h=0.8, label="hot field, no fluid exchange", kind="hot")],
                  wells=[dict(x=0.45, z=3.0, kind="coaxial", dir="loop", label="coaxial, VIT")],
                  units=[dict(x=0.74, label="binary plant", kind="orc")],
                  flows=[dict(a="w0", b="u0", label="56–60 °C", hot=True), dict(a="u0", b="w0", label="40 °C", hot=False, back=True)])),
    8: dict(
        short="Multilateral deep closed loop",
        sheet="Two vertical wells joined by many long horizontal laterals (Eavor-Loop), water circulated in a closed loop by thermosiphon; "
              "heat by conduction only. The source table carries 8.2 MWe (company) against 2.2 MWe (independent) per loop, 64 MWth.",
        example="Eavor-Loop 1 type: two 4.5 km verticals, 12 laterals of 5 km (69 km of hole in all) in 55 °C/km rock (263 °C at depth), "
                "8½ in open-hole laterals 80 m apart, 100 kg/s injected at 40 °C against 2 MPa.",
        model="Injector → lateral stack → producer as a chain of segments with transient hand-over; laterals in an adiabatic unit cell "
              "80 m wide (and, as a bound, isolated laterals in an infinite rock); production wellhead pressure from the march.",
        look="The laterals start to interfere after a few years; the 80 m spacing costs about 5 K by year 30.",
        plant="binary (GEOPHIRES practical)",
        spec=dict(depth=5.0, ticks=[0, 2, 4], isotherms=[(4.5, "263 °C")],
                  wells=[dict(x=0.22, z=4.5, kind="single", dir="down", label="injector", lateral=dict(x2=0.78, n=12)),
                         dict(x=0.78, z=4.5, kind="single", dir="up", label="producer")],
                  units=[dict(x=0.5, label="binary plant", kind="orc")],
                  flows=[dict(a="w1", b="u0", label="111–136 °C", hot=True), dict(a="u0", b="w0", label="40 °C", hot=False)],
                  notes=[dict(x=0.5, z=4.1, text="12 laterals × 5 km, 80 m apart")])),
    9: dict(
        short="Doublet with horizontal closed-loop connections",
        sheet="A doublet whose wells are connected by one or more horizontal closed-loop sections (advanced closed loop, Chevron / Baseload "
              "type); no public parameters. The source table lists capacity as undisclosed.",
        example="Assumed: 3.0 km verticals, two 3 km laterals 120 m apart, 40 °C/km (135 °C at depth), 35 kg/s at 40 °C.",
        model="Same chain as type 8 with two laterals in a 120 m unit cell.",
        look="A short loop in warm rock: 1.8 MWth and 52 °C water after 30 years, useful for heat only.",
        plant="binary (GEOPHIRES practical)",
        spec=dict(depth=3.4, ticks=[0, 1, 2, 3], isotherms=[(3.0, "135 °C")],
                  wells=[dict(x=0.25, z=3.0, kind="single", dir="down", label="injector", lateral=dict(x2=0.75, n=2)),
                         dict(x=0.75, z=3.0, kind="single", dir="up", label="producer")],
                  units=[dict(x=0.5, label="heat / binary", kind="hx")],
                  flows=[dict(a="w1", b="u0", label="52–57 °C", hot=True), dict(a="u0", b="w0", label="40 °C", hot=False)],
                  notes=[dict(x=0.5, z=2.6, text="2 laterals × 3 km (assumed)")])),
    10: dict(
        short="Doublet EGS in hot dry rock",
        sheet="Engineered geothermal systems in crystalline rock: an injector and a producer connected through stimulated fractures "
              "(Soultz, United Downs). The source table quotes 3–5 MWe and about 12 MWth per well.",
        example="5.0 km doublet in 200 °C granite: 20 stimulated fractures 50 m apart, 600 m between the wells, 400 m fracture breadth, "
                "30 kg/s, 60 °C injection, impedance 0.3 MPa per L/s, production pump at 500 m, 2 MPa wellhead.",
        model="Injector and producer in single-channel modes with the Gringarten (1975) parallel-fracture reservoir between them (as in "
              "GEOPHIRES), Stehfest inversion each time step; impedance charged as circulation pump work.",
        look="With 20 well-spaced fractures the thermal front needs about 140 years to cross 600 m, so the reservoir holds 200 °C for the "
             "whole 30 years and the produced water warms slightly as the casing heats up.",
        plant="binary (GEOPHIRES practical)",
        spec=dict(depth=5.6, ticks=[0, 2, 4], isotherms=[(5.0, "200 °C")],
                  wells=[dict(x=0.28, z=5.0, kind="single", dir="down", label="injector"),
                         dict(x=0.72, z=5.0, kind="single", dir="up", label="producer", pump=0.5)],
                  fractures=dict(x1=0.28, x2=0.72, z=5.0, n=20, h=0.45),
                  units=[dict(x=0.5, label="binary plant", kind="orc")],
                  flows=[dict(a="w1", b="u0", label="182–185 °C", hot=True), dict(a="u0", b="w0", label="60 °C", hot=False)],
                  notes=[dict(x=0.5, z=4.35, text="20 fractures, 600 m apart")])),
    11: dict(
        short="Hybrid single-well EGS",
        sheet="A single deep well with injection into an upper fracture zone and production from a lower one through a short fractured "
              "path, heat exchanged at surface (Eden Project type). The source table lists capacity as small.",
        example="4.8 km coaxial well in 35 °C/km rock (183 °C at TD): injection down the annulus, 150 m path through 3 fractures 30 m apart, "
                "return up 3½ × 4½ in insulated tubing, 10 kg/s at 40 °C.",
        model="Coaxial well with a <code>bottom_element</code>: a parallel-fracture reservoir in series at TD, its outlet enthalpy handed "
              "to the riser; the fluid also exchanges heat with the rock along the whole well.",
        look="A small fractured volume: the produced water falls from 133 °C in year 1 to 57 °C after 30 years.",
        plant="binary (GEOPHIRES practical)",
        spec=dict(depth=5.3, ticks=[0, 2, 4], isotherms=[(4.8, "183 °C")],
                  wells=[dict(x=0.45, z=4.8, kind="coaxial", dir="loop", label="single well", open_bottom=True)],
                  fractures=dict(x1=0.33, x2=0.57, z=4.75, n=3, h=0.35),
                  units=[dict(x=0.74, label="heat exchanger", kind="hx")],
                  flows=[dict(a="w0", b="u0", label="57–133 °C", hot=True), dict(a="u0", b="w0", label="40 °C", hot=False, back=True)],
                  notes=[dict(x=0.45, z=4.2, text="3 fractures, 150 m path")])),
    12: dict(
        short="Single-well cyclic fracture loop",
        sheet="Pressure geothermal (Sage HeatRoot type): one well with one large propped fracture, operated in daily cycles of injection, "
              "soak and production. The source table gives no capacity.",
        example="3.5 km well in 150 °C rock, one fracture 300 m × 100 m × 8 mm; each day 200 m³ injected at 40 °C over 4 h, 12 h soak, "
                "produced over 8 h; 365 cycles, with the 3.5 km production leg modelled.",
        model="<code>CyclicFracture</code>: plug flow in the fracture (last in, first out) with implicit 1-D conduction into the rock on each "
              "cell, checked against Lauwerier; the daily-mean produced temperature fed to a production-well model.",
        look="The first cycles return 84 °C water; the fracture's neighbourhood cools and the output settles near 72 °C.",
        plant="binary on the daily-average stream",
        spec=dict(depth=3.9, ticks=[0, 1, 2, 3], isotherms=[(3.5, "150 °C")],
                  wells=[dict(x=0.45, z=3.5, kind="single", dir="cycle", label="inject / soak / produce")],
                  fractures=dict(x1=0.3, x2=0.6, z=3.5, n=1, h=0.5),
                  units=[dict(x=0.76, label="binary (daily batches)", kind="orc")],
                  flows=[dict(a="w0", b="u0", label="72–84 °C, 8 h a day", hot=True), dict(a="u0", b="w0", label="40 °C, 4 h a day", hot=False, back=True)],
                  notes=[dict(x=0.45, z=2.9, text="one fracture 300 × 100 m")])),
    13: dict(
        short="Multi-well cyclic EGS with storage",
        sheet="Several cyclic wells phased for continuous output, the fractures also used as thermal storage charged with surplus "
              "electricity (Sage HeatCycle / Battery+ type). The source table quotes up to 50 MWe per 18–20 wells.",
        example="Twenty type-12 wells phased through the day; storage test: a cycled fracture charged with 200 m³ of 150 °C water, 24 h soak, "
                "discharge, compared with the uncharged cycle.",
        model="Type 12 scaled by 20 for the continuous field; the storage round trip from the same cyclic-fracture model (energy recovered "
              "above the baseline cycle, then converted at the ORC efficiency for the discharge temperature).",
        look="Thermal round trip of the stored heat and the small electric round trip once it is converted back.",
        plant="binary; storage round trip",
        spec=dict(depth=3.9, ticks=[0, 1, 2, 3], isotherms=[(3.5, "150 °C")],
                  wells=[dict(x=0.2, z=3.5, kind="single", dir="cycle", label=""), dict(x=0.38, z=3.5, kind="single", dir="cycle", label=""),
                         dict(x=0.56, z=3.5, kind="single", dir="cycle", label="20 wells phased")],
                  fractures=dict(at=[0.2, 0.38, 0.56], z=3.5, n=3, h=0.5),
                  units=[dict(x=0.84, label="binary + heaters", kind="heater")],
                  flows=[dict(a="w2", b="u0", label="72 °C continuous", hot=True), dict(a="u0", b="w2", label="40 °C / 150 °C charge", hot=False, back=True)])),
    14: dict(
        short="Horizontal doublet EGS",
        sheet="Horizontal injector and producer with many plug-and-perf fracture stages between them (Fervo Project Red, Cape Station). "
              "The source table quotes the 3.5 MWe pilot at 63 L/s.",
        example="Project Red type: 2.3 km TVD plus 1 km laterals 100 m apart in 190 °C rock, 16 stages with 4 clusters each (64 fractures "
                "15.6 m apart, 250 m breadth), 58 kg/s (63 L/s), 60 °C injection, ESP at 1000 m, 2 MPa wellhead, impedance 0.15 MPa per L/s.",
        model="Deviated injector and producer (survey with the horizontal section) with the parallel-fracture reservoir between them; "
              "the ESP sits at 1000 m because the stream would flash below a 500 m pump.",
        look="Year 1 reproduces the pilot's 3.5 MWe; the rock between laterals 100 m apart is swept out in about ten years.",
        plant="binary (GEOPHIRES practical)",
        spec=dict(depth=2.7, ticks=[0, 1, 2], isotherms=[(2.3, "190 °C")],
                  wells=[dict(x=0.2, z=2.3, kind="single", dir="down", label="injector", lateral=dict(x2=0.6, n=1, dy=-0.1)),
                         dict(x=0.32, z=2.3, kind="single", dir="up", label="producer", pump=1.0, lateral=dict(x2=0.72, n=1, dy=0.1))],
                  fractures=dict(x1=0.36, x2=0.7, z=2.3, n=16, h=0.25, horizontal_pair=True),
                  units=[dict(x=0.82, label="binary plant", kind="orc")],
                  flows=[dict(a="w1", b="u0", label="184 °C → 64 °C", hot=True), dict(a="u0", b="w0", label="60 °C", hot=False)],
                  notes=[dict(x=0.53, z=1.75, text="laterals 100 m apart, 64 fractures")])),
    15: dict(
        short="Co-produced water from oil and gas wells",
        sheet="Hot water that producing oil and gas wells already lift, separated at surface and passed through a binary unit before "
              "disposal. The source table quotes 20–350 kWe per field.",
        example="Eight 2.0 km wells producing 15 kg/s of liquid each from a 95 °C reservoir up 2⅞ in tubing with a static water annulus "
                "inside 7 in casing; the field's 120 kg/s to an ORC returning at 60 °C.",
        model="One production well in tubing mode (static annulus as a conducting layer), constant 95 °C reservoir; the lifting power is part "
              "of oil production and not charged.",
        look="The small tubing loses little heat: 93 °C at surface, 0.72 MWe for the eight wells.",
        plant="binary (GEOPHIRES practical)",
        spec=dict(depth=2.3, ticks=[0, 1, 2], isotherms=[(2.0, "95 °C")],
                  bands=[dict(z=1.9, h=0.2, label="oil reservoir", kind="aquifer")],
                  wells=[dict(x=0.22, z=2.0, kind="single", dir="up", label="8 oil wells", tubing=True, label_side="left"), dict(x=0.36, z=2.0, kind="single", dir="up", label="", tubing=True),
                         dict(x=0.5, z=2.0, kind="single", dir="up", label="", tubing=True), dict(x=0.86, z=2.0, kind="single", dir="down", label="disposal")],
                  units=[dict(x=0.68, label="separator → ORC", kind="orc")],
                  flows=[dict(a="w2", b="u0", label="93 °C water", hot=True), dict(a="u0", b="w3", label="60 °C", hot=False)])),
    16: dict(
        short="Closed loop in an abandoned well",
        sheet="A single abandoned hydrocarbon well converted to a coaxial closed loop (CeraPhi, wells2watts type) for heat. The source table "
              "quotes under 1 MWth (a 400 kW trial).",
        example="3 km gas well, 7 in casing, 110 °C at TD, 3½ × 4½ in insulated tubing, 5 kg/s at 30 °C inlet, 1 MPa.",
        model="The coaxial solver, 30 years.",
        look="39–42 °C water and 0.19 MWth: a heat source for a building, not a power source.",
        plant="heat only",
        spec=dict(depth=3.4, ticks=[0, 1, 2, 3], isotherms=[(3.0, "110 °C")],
                  wells=[dict(x=0.45, z=3.0, kind="coaxial", dir="loop", label="abandoned well, VIT")],
                  units=[dict(x=0.74, label="building heat", kind="heat")],
                  flows=[dict(a="w0", b="u0", label="39–42 °C", hot=True), dict(a="u0", b="w0", label="30 °C", hot=False, back=True)])),
    17: dict(
        short="Heat recovery after thermal EOR",
        sheet="Recovering the heat left in a steam-flooded (SAGD) reservoir once oil production ends, by circulating water through the hot "
              "chamber with high-temperature ESPs. The source table lists it as to be proven.",
        example="A SAGD chamber 1000 × 100 × 30 m at 300 m depth, 200 °C, swept by a doublet at 20 kg/s with 50 °C return, 15 years.",
        model="<code>HeatStore</code>: lumped chamber with Duhamel conductive losses to 12 °C surroundings; the 300 m wells lose little and "
              "are not modelled.",
        look="The chamber holds a few years of heat: 112 °C in year 1, below the return temperature by year 7.",
        plant="binary (GEOPHIRES practical)",
        spec=dict(depth=0.36, ticks=[0, 0.1, 0.2, 0.3], isotherms=[(0.3, "200 °C → ambient")],
                  bands=[dict(z=0.285, h=0.03, label="steam chamber 1000 × 100 × 30 m", kind="store")],
                  wells=[dict(x=0.3, z=0.3, kind="single", dir="up", label="producer (HT ESP)", pump=0.15),
                         dict(x=0.7, z=0.3, kind="single", dir="down", label="injector")],
                  units=[dict(x=0.5, label="binary / heat", kind="orc")],
                  flows=[dict(a="w0", b="u0", label="112 → 58 °C", hot=True), dict(a="u0", b="w1", label="50 °C", hot=False)])),
    18: dict(
        short="Supercritical fluid below a volcanic field",
        sheet="Wells drilled below conventional reservoirs into supercritical conditions (IDDP-1, IDDP-2, Japan Beyond-Brittle); the fluid "
              "flashes to superheated steam on the way up. The source table expects over 30 MWe and over 75 MWth per well.",
        example="IDDP-2 type: 4.5 km in a 100 °C/km field, 450 °C / 30 MPa feed (productivity index 10 kg/s/MPa), 20 kg/s to a 10 MPa "
                "wellhead; 9⅝ in casing to 3 km, 7 in below.",
        model="Single production well fed by enthalpy (IAPWS-95 to 1000 MPa / 2000 K tables), supercritical-to-steam expansion in the riser "
              "with the enthalpy-step and relaxation safeguards; condensing turbine on the superheated stream.",
        look="The well delivers steam at the 10 MPa saturation temperature (311 °C) with a small superheat: 13.4 MWe from 20 kg/s.",
        plant="superheated-steam turbine",
        spec=dict(depth=5.0, ticks=[0, 2, 4], isotherms=[(4.5, "450 °C"), (2.0, "215 °C")],
                  bands=[dict(z=4.2, h=0.8, label="supercritical zone, 30 MPa", kind="hot"), dict(z=1.6, h=0.8, label="conventional field", kind="steam")],
                  wells=[dict(x=0.42, z=4.5, kind="single", dir="up", label="IDDP type", steam=True)],
                  units=[dict(x=0.72, label="steam turbine", kind="steam")],
                  flows=[dict(a="w0", b="u0", label="311 °C steam, 10 MPa", hot=True)])),
    19: dict(
        short="Superhot-rock EGS doublet",
        sheet="EGS in rock above 375 °C (Newberry, AltaRock / Mazama, Quaise's target): fractured superhot rock produces steam. The source "
              "table quotes up to 50 MWe and 125 MWth per well, marked as uncertain.",
        example="3.5 km doublet in 375 °C rock, 20 fractures 40 m apart, 300 m path, 400 m breadth, 40 kg/s, 60 °C injection at 5 MPa, "
                "production up 7 in casing to a 5 MPa wellhead.",
        model="Parallel-fracture reservoir between single-channel wells; the producer flashes in the riser to a two-phase wellhead stream "
              "(steam fraction 0.35 in year 1); double flash at the wellhead pressure.",
        look="The fractured block cools from 375 to 285 °C in 30 years; the wellhead stays at saturation for 5 MPa while the steam "
             "fraction falls from 0.35 to 0.03.",
        plant="double flash",
        spec=dict(depth=4.0, ticks=[0, 1, 2, 3], isotherms=[(3.5, "375 °C")],
                  wells=[dict(x=0.28, z=3.5, kind="single", dir="down", label="injector"),
                         dict(x=0.72, z=3.5, kind="single", dir="up", label="producer", flash=1.5)],
                  fractures=dict(x1=0.28, x2=0.72, z=3.5, n=20, h=0.45),
                  units=[dict(x=0.5, label="flash plant", kind="flash")],
                  flows=[dict(a="w1", b="u0", label="264 °C, x 0.35 → 0.03", hot=True), dict(a="u0", b="w0", label="60 °C", hot=False)],
                  notes=[dict(x=0.5, z=2.9, text="20 fractures, 300 m apart")])),
    20: dict(
        short="Ultradeep millimetre-wave drilled loop",
        sheet="Boreholes to 10–20 km drilled with millimetre-wave energy into 500 °C rock (Quaise type), vitrified walls, operated here as the "
              "source table classes it: an ultradeep closed loop. The source table quotes 100 MWe per plant.",
        example="20 km coaxial loop, 500 °C at TD (24.25 °C/km), 0.30 m glass-lined bore (2 cm glass, k 1.4 W/m/K), 5½ × 7 in insulated "
                "tubing, 20 and 40 kg/s injected at 40 °C against 10 MPa.",
        model="The coaxial solver with the glass liner as a completion layer, 100 m cells, 30 years; binary, flash and steam options compared.",
        look="Conduction-starved like every closed loop: 171 °C falling to 146 °C at 20 kg/s, and less power at 40 kg/s than at 20.",
        plant="binary (best of the options)",
        spec=dict(depth=22.0, ticks=[0, 5, 10, 15, 20], isotherms=[(20, "500 °C"), (10, "258 °C")],
                  bands=[dict(z=18.5, h=3.0, label="hot basement, no fluid exchange", kind="hot")],
                  wells=[dict(x=0.45, z=20.0, kind="coaxial", dir="loop", label="20 km coaxial, glass-lined")],
                  units=[dict(x=0.74, label="binary plant", kind="orc")],
                  flows=[dict(a="w0", b="u0", label="146–171 °C", hot=True), dict(a="u0", b="w0", label="40 °C", hot=False, back=True)])),
}

SHEET_QUOTE = {   # capacity quoted in the source table (electrical, thermal, basis)
    1: ("<2 MWe", "<5 MWth", "per well"), 2: ("2–6 MWe (<7)", "<15 MWth", "per well"), 3: ("3–25 MWe (up to 30–50)", "12.5 MWth typical; 70–125 MWth", "per well"),
    3.5: ("3–25 MWe", "–", "per well"), 4: ("355 kWe per module; 1–2 MWe stacked", "cooled brine to district heating", "per module"),
    5: ("—", "<10 kW", "per unit"), 6: ("—", "<10 kW to <5 MWth", "per installation"), 7: ("2–9 MWe claimed; 2 MWe retrofit", "heat for buildings / data centres", "per well"),
    8: ("2.2 MWe (independent) vs 8.2 MWe (company)", "64 MWth (company)", "per loop"), 9: ("undisclosed", "undisclosed", "–"),
    10: ("3–5 MWe", "~12 MWth", "per well"), 11: ("small", "small", "per well"), 12: ("not stated", "not stated", "per well"),
    13: ("up to 50 MWe", "—", "per system (18–20 wells)"), 14: ("3.5 MWe pilot (63 L/s)", "heat (300,000 homes planned)", "per pilot"),
    15: ("20–350 kWe", "heating / cooling", "per field"), 16: ("—", "<1 MWth (~400 kW trial)", "per well"), 17: ("to be proven", "to be proven", "–"),
    18: (">30 MWe expected", ">75 MWth expected", "per well"), 19: ("up to 50 MWe (??)", "up to 125 MWth (??)", "per well"), 20: ("100 MWe", "—", "per plant"),
}

# numeric range of the quoted electrical capacity for the comparison chart (MWe, low, high), None where no figure
SHEET_MWE = {1: (0.0, 2.0), 2: (2.0, 7.0), 3: (3.0, 25.0), 3.5: (3.0, 25.0), 4: (0.355, 2.0), 5: None, 6: None, 7: (2.0, 9.0), 8: (2.2, 8.2),
             9: None, 10: (3.0, 5.0), 11: None, 12: None, 13: (50.0, 50.0), 14: (3.5, 3.5), 15: (0.02, 0.35), 16: None, 17: None,
             18: (30.0, 30.0), 19: (50.0, 50.0), 20: (100.0, 100.0)}

TINY = {1: "HSA doublet", 2: "Pumped MT well", 3: "HT flash well", 3.5: "Dry-steam well", 4: "Waste-heat ORC", 5: "GSHP U-tube", 6: "400 m coaxial BHE",
        7: "GreenLoop retrofit", 8: "Eavor-type loop", 9: "Closed-loop doublet", 10: "EGS doublet", 11: "Single-well EGS", 12: "Cyclic fracture",
        13: "20 cyclic wells", 14: "Fervo-type doublet", 15: "Co-produced water", 16: "Abandoned-well loop", 17: "SAGD heat store",
        18: "Supercritical well", 19: "Superhot EGS", 20: "20 km loop"}
