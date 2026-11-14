"""
12b_presentation.py

PowerPoint for the skipper: what the sonar showed (ST90 and SU90, SCH302).
Uses the figures and videos of 12_presentation_figures.R plus two existing
animations. Run from the project root after that script:
    python 12b_presentation.py
Output: output_presentation/SCH302_sonar_results_for_skipper.pptx
"""
import os
from PIL import Image
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.chart.data import CategoryChartData
from pptx.enum.chart import XL_CHART_TYPE, XL_LABEL_POSITION

OUT = "output_presentation"
VID = {"follow": "output_ST90/animation/vessel_schools_follow_skip1000m_2026-07-02.mp4",
       "origin": "output_ST90/animation/schools_origin.mp4",
       "drift": os.path.join(OUT, "H_herring_drift.mp4")}

NAVY, ORANGE, TEAL = RGBColor(0x16, 0x32, 0x4F), RGBColor(0xE8, 0x83, 0x3A), RGBColor(0x2A, 0x9D, 0xB5)
INK, MUTED, TINT, WHITE = RGBColor(0x1E, 0x29, 0x33), RGBColor(0x5B, 0x67, 0x73), RGBColor(0xEE, 0xF5, 0xF8), RGBColor(0xFF, 0xFF, 0xFF)
HEAD, BODY = "Cambria", "Calibri"

prs = Presentation()
prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
BLANK = prs.slide_layouts[6]
SW, SH = 13.333, 7.5


def text(slide, x, y, w, h, runs, size=18, color=INK, font=BODY, bold=False, align=PP_ALIGN.LEFT,
         anchor=MSO_ANCHOR.TOP, space_after=6, name=None):
    """runs: str, or list of paragraphs; a paragraph is a str or a list of (text, {bold, color, size}) runs."""
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    if name:
        tb.name = name
    tf = tb.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    paras = runs if isinstance(runs, list) else [runs]
    for i, para in enumerate(paras):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.space_after = Pt(space_after)
        for t, o in ([(para, {})] if isinstance(para, str) else para):
            r = p.add_run()
            r.text = t
            r.font.name = o.get("font", font)
            r.font.size = Pt(o.get("size", size))
            r.font.bold = o.get("bold", bold)
            r.font.italic = o.get("italic", False)
            r.font.color.rgb = o.get("color", color)
    return tb


def title(slide, t, dark=False):
    text(slide, 0.6, 0.4, SW - 1.2, 0.9, t, size=34, font=HEAD, bold=True, color=WHITE if dark else NAVY,
         anchor=MSO_ANCHOR.MIDDLE, name="Title")


def fit_image(slide, path, x, y, w, h, name=None):
    iw, ih = Image.open(path).size
    s = min(w / iw, h / ih)
    pw, ph = iw * s, ih * s
    pic = slide.shapes.add_picture(path, Inches(x + (w - pw) / 2), Inches(y + (h - ph) / 2), Inches(pw), Inches(ph))
    if name:
        pic.name = name
    return pic


def movie(slide, path, poster, x, y, w, h):
    iw, ih = Image.open(poster).size
    s = min(w / iw, h / ih)
    pw, ph = iw * s, ih * s
    return slide.shapes.add_movie(path, Inches(x + (w - pw) / 2), Inches(y + (h - ph) / 2), Inches(pw), Inches(ph),
                                  poster_frame_image=poster, mime_type="video/mp4")


def shape(slide, kind, x, y, w, h, fill, line=None, name=None):
    s = slide.shapes.add_shape(kind, Inches(x), Inches(y), Inches(w), Inches(h))
    s.fill.solid()
    s.fill.fore_color.rgb = fill
    if line is None:
        s.line.fill.background()
    else:
        s.line.color.rgb = line
    s.shadow.inherit = False
    if name:
        s.name = name
    return s


def badge(slide, x, y, label, d=0.55, fill=ORANGE):
    c = shape(slide, MSO_SHAPE.OVAL, x, y, d, d, fill)
    tf = c.text_frame
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.CENTER
    r = p.add_run()
    r.text = label
    r.font.size, r.font.bold, r.font.name, r.font.color.rgb = Pt(20), True, BODY, WHITE
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    return c


def card(slide, x, y, w, h, head, body, accent=ORANGE, big=None):
    shape(slide, MSO_SHAPE.ROUNDED_RECTANGLE, x, y, w, h, TINT).adjustments[0] = 0.08
    yy = y + 0.2
    if big:
        text(slide, x + 0.25, yy, w - 0.5, 0.75, big, size=36, bold=True, color=accent, font=HEAD)
        yy += 0.8
    text(slide, x + 0.25, yy, w - 0.5, 0.45, head, size=18, bold=True, color=NAVY)
    text(slide, x + 0.25, yy + 0.45, w - 0.5, h - (yy - y) - 0.55, body, size=15, color=MUTED)


def notes(slide, t):
    slide.notes_slide.notes_text_frame.text = t


def caption(slide, t, y=6.85):
    text(slide, 0.6, y, SW - 1.2, 0.4, t, size=12, color=MUTED)


def dark_bg(slide):
    bg = slide.background.fill
    bg.solid()
    bg.fore_color.rgb = NAVY


# ---------------------------------------------------------------------------------------------- 1 title
s = prs.slides.add_slide(BLANK)
dark_bg(s)
text(s, 0.8, 2.0, 6.8, 1.6, "What the sonar saw", size=54, font=HEAD, bold=True, color=WHITE)
text(s, 0.8, 3.3, 6.6, 1.4, "Herring and mackerel schools around SCH302, from the ST90 and SU90 sonars", size=24, color=RGBColor(0xCF, 0xDD, 0xE8))
text(s, 0.8, 5.6, 7.0, 0.5, "October 2025 and June–July 2026", size=18, color=ORANGE, bold=True)
fit_image(s, os.path.join(OUT, "poster_drift.png"), 7.9, 1.0, 5.0, 5.5)
notes(s, "Thank you for letting us use the sonar data. Today: what we could learn from the schools the sonar picked up, "
         "how they move, and whether that says anything about catches. Questions are welcome at any point.")

# ---------------------------------------------------------------------------------------------- 2 where
s = prs.slides.add_slide(BLANK)
title(s, "Three sets of sonar recordings east of Shetland")
fit_image(s, os.path.join(OUT, "A_where.png"), 0.5, 1.4, 6.6, 5.5)
card(s, 7.5, 1.45, 5.3, 1.55, "October 2025 · SU90 · mackerel", "23–29 October, 13 hauls, 1,175 schools detected", accent=TEAL)
card(s, 7.5, 3.2, 5.3, 1.55, "June–July 2026 · ST90 · herring", "28 June – 5 July, 22 hauls, about 6,000 schools detected")
card(s, 7.5, 4.95, 5.3, 1.55, "June–July 2026 · SU90 · herring", "Same weeks, both sonars running side by side: a rare chance to compare them")
notes(s, "Blue dots: where the SU90 picked up mackerel schools in October 2025. Orange: where the ST90 picked up herring schools "
         "in June–July 2026. In 2026 both sonars were running at the same time, so we can compare them.")

# ---------------------------------------------------------------------------------------------- 3 how
s = prs.slides.add_slide(BLANK)
title(s, "From echoes to fish movement in four steps")
steps = [("1", "Detect", "The sonar software (LSSS Profos) picks out schools on every ping and gives their position"),
         ("2", "Clean", "Remove echoes that are not fish, for example echoes that travel with the ship"),
         ("3", "Smooth", "Positions jump a few metres from ping to ping; smoothing gives each school's real track and speed"),
         ("4", "Compare", "Put the tracks next to the catches, the tide and the current model")]
for i, (n, h, b) in enumerate(steps):
    x = 0.6 + i * 3.1
    shape(s, MSO_SHAPE.ROUNDED_RECTANGLE, x, 2.0, 2.8, 3.9, TINT).adjustments[0] = 0.06
    badge(s, x + 0.25, 2.25, n, d=0.75)
    text(s, x + 0.25, 3.2, 2.35, 0.5, h, size=24, bold=True, color=NAVY, font=HEAD)
    text(s, x + 0.25, 3.8, 2.35, 2.0, b, size=16, color=MUTED)
    if i < 3:
        a = shape(s, MSO_SHAPE.RIGHT_ARROW, x + 2.83, 3.7, 0.25, 0.4, ORANGE)
text(s, 0.6, 6.25, 12, 0.6, "Every result in this talk comes from thousands of schools, not from one or two.", size=18, color=NAVY, bold=True)
notes(s, "Keep this short. The point: the sonar gives a position for each school on each ping. We clean out false echoes, "
         "smooth the jumpy positions into a track, and then compare with catch, tide and currents.")

# ---------------------------------------------------------------------------------------------- 4 cleaning
s = prs.slides.add_slide(BLANK)
title(s, "Some 'schools' were really the ship's own echoes")
fit_image(s, os.path.join(OUT, "B_cleaning.png"), 0.4, 1.4, 8.4, 5.3)
card(s, 9.1, 1.6, 3.7, 2.2, "Real schools", "Mostly under 2 knots", big="~1 kn")
card(s, 9.1, 4.0, 3.7, 2.6, "Removed", "Echoes moving at the ship's steaming speed and heading, plus other unrealistic tracks",
     big="1 in 9")
notes(s, "When we followed each 'school', a group moved at 12–14 knots, in the same direction as the ship. Fish do not do that: "
         "these are echoes tied to the vessel (wake, bubbles, the hull). We removed them, about 1 in 9 tracks on the ST90. "
         "The real schools mostly move at about 1 knot.")

# ---------------------------------------------------------------------------------------------- 5 smoothing
s = prs.slides.add_slide(BLANK)
title(s, "Smoothing turns jumpy dots into a real track")
fit_image(s, os.path.join(OUT, "C_smoothing.png"), 0.4, 1.4, 7.6, 5.4)
text(s, 8.4, 1.8, 4.4, 4.5, [
    [("Grey dots: ", {"bold": True, "color": NAVY}), ("where the sonar put the school on each ping. They jump 3–5 m every time.", {})],
    [("Orange line: ", {"bold": True, "color": ORANGE}), ("the smoothed track of the same school.", {})],
    [("Why it matters: ", {"bold": True, "color": NAVY}), ("measured from the raw dots, every school would seem to dart about at several knots. "
                                                         "The smoothed track gives its real speed, here 1.2 knots.", {})]],
     size=18, space_after=14)
notes(s, "One herring school followed for 7 minutes. The dots jitter; the line is what the school really did.")

# ---------------------------------------------------------------------------------------------- 6 video: around the vessel
s = prs.slides.add_slide(BLANK)
title(s, "The schools around the vessel during a tow")
movie(s, VID["follow"], os.path.join(OUT, "poster_follow.png"), 0.5, 1.35, 6.6, 5.6)
text(s, 7.6, 1.8, 5.2, 4.8, [
    [("2 July, ST90, haul 53. ", {"bold": True, "color": NAVY}), ("Click the picture to play.", {"color": MUTED})],
    "The red boat is the vessel, always in the middle, north up. Rings are 200, 400 and 600 m.",
    "Dots are schools; the tails show where they are swimming. School movement is shown 5 times bigger and 10 times slower, so it can be seen.",
    "Stretches with no school within 1 km of the vessel are skipped."],
     size=18, space_after=14)
notes(s, "Let the video run. Point out schools crossing ahead of the vessel, and that they move in different directions.")

# ---------------------------------------------------------------------------------------------- 7 video: every school its own way
s = prs.slides.add_slide(BLANK)
title(s, "Each school goes its own way")
movie(s, VID["origin"], os.path.join(OUT, "poster_origin.png"), 0.5, 1.35, 6.6, 5.6)
text(s, 7.6, 1.8, 5.2, 4.8, [
    "Here every school is moved to the same starting point and set off at the same moment.",
    "They spread out in all directions, like fireworks.",
    [("So: ", {"bold": True, "color": NAVY}), ("watching a few schools tells you little about where the herring as a whole are going. "
                                              "For that we need the average of thousands of schools, and we need to remove the tide.", {})]],
     size=18, space_after=14)
notes(s, "This sets up the next slides: individual schools wander; any common direction is a small push on top of that.")

# ---------------------------------------------------------------------------------------------- 8 tide
s = prs.slides.add_slide(BLANK)
title(s, "The tide pushes the schools back and forth")
fit_image(s, os.path.join(OUT, "E_tide.png"), 0.4, 1.35, 8.6, 5.3)
card(s, 9.3, 1.6, 3.5, 2.2, "Tidal push", "On the schools, each way", big="0.3 kn")
card(s, 9.3, 4.0, 3.5, 2.6, "Turns every ~6 hours", "Over a day it cancels out, so we remove it to see what is left over", accent=TEAL,
     big="~6 h")
caption(s, "Dots: average north–south speed of the herring schools each hour (ST90, 1–2 July). Blue line: the tide fitted to all 8 days.")
notes(s, "The schools are carried north, then south, with the tide. This is bigger than any steady movement, "
         "so you only see the steady movement after taking the tide out.")

# ---------------------------------------------------------------------------------------------- 9 drift video
s = prs.slides.add_slide(BLANK)
title(s, "With the tide removed, herring drifted west-south-west")
movie(s, VID["drift"], os.path.join(OUT, "poster_drift.png"), 0.4, 1.3, 7.6, 5.35)
card(s, 8.4, 1.5, 4.4, 2.3, "Average daylight drift", "Towards west-south-west (237°)", big="≈ 5 nm a day")
text(s, 8.4, 4.1, 4.4, 2.8, [
    "The orange line adds up the average daytime movement of all herring schools, hour after hour (tide and ship effect removed).",
    "The same direction comes out whichever day we leave out."],
     size=16, color=MUTED, space_after=10)
caption(s, "Only daylight hours are joined together, so the line shows the trend, not the path of real fish.")
notes(s, "About 0.2 knots, or roughly 5 nautical miles a day, towards Fair Isle and Orkney. Small compared with the tide, "
         "but the same every day. The line joins only the daylight hours, so it is a picture of the trend.")

# ---------------------------------------------------------------------------------------------- 10 current or swimming
s = prs.slides.add_slide(BLANK)
title(s, "Was it the current, or were the fish swimming?")
fit_image(s, os.path.join(OUT, "F_current_or_swimming.png"), 0.3, 1.25, 8.6, 5.5)
text(s, 9.1, 1.5, 3.8, 5.3, [
    [("Herring 2026: ", {"bold": True, "color": ORANGE}), ("the current (Copernicus ocean model) was weak and ran south-east. "
                                                          "So the herring were swimming west-south-west themselves, about 5–6 nm a day.", {})],
    [("Mackerel 2025: ", {"bold": True, "color": TEAL}), ("they also went south, but the current was running the same way. "
                                                         "They were mostly carried.", {})],
    [("Herring direction ", {"bold": True, "color": NAVY}), ("points towards the Orkney–Shetland spawning grounds (August–September). "
                                                            "It fits a pre-spawning movement, but does not prove it.", {})]],
     size=16, space_after=14)
notes(s, "Arrows: dark = movement over the seabed, blue = current, orange = swimming through the water (dark minus blue). "
         "The comparison with the current is what tells the two years apart.")

# ---------------------------------------------------------------------------------------------- 11 catch
s = prs.slides.add_slide(BLANK)
title(s, "More schools on the sonar did not mean a bigger catch")
fit_image(s, os.path.join(OUT, "G_catch.png"), 0.4, 1.35, 7.4, 5.3)
reasons = [("1", "The sonar sees only part of it", "Ahead of the vessel, out to about 600 m"),
           ("2", "Catches are capped", "By your decisions and the vessel's capacity: half of the hauls were 110–160 t"),
           ("3", "Counts describe the area", "Schools within 5 km, not the fish that went into the net")]
for i, (n, h, b) in enumerate(reasons):
    y = 1.6 + i * 1.65
    badge(s, 8.2, y, n)
    text(s, 9.0, y - 0.05, 3.9, 0.45, h, size=18, bold=True, color=NAVY)
    text(s, 9.0, y + 0.42, 3.9, 0.9, b, size=15, color=MUTED)
caption(s, "Each dot is one herring haul in 2026 (number = haul). The same holds for the SU90 and for the mackerel in 2025.")
notes(s, "No pattern: hauls with many schools around them were not bigger. Your experience is probably the better guide here; "
         "the sonar is more useful for seeing how fish behave.")

# ---------------------------------------------------------------------------------------------- 12 depth
s = prs.slides.add_slide(BLANK)
title(s, "The depth in the sonar export is not the school's depth")
fit_image(s, os.path.join(OUT, "D_beams.png"), 0.4, 1.35, 8.4, 5.2)
card(s, 9.1, 1.6, 3.7, 2.3, "Same school, two depths", "ST90 says 28 m, SU90 says 63 m", big="28 vs 63 m")
text(s, 9.1, 4.2, 3.7, 2.4, "The software gives the depth of the centre of the beam at that distance. "
                            "The real depth needs the echosounder.", size=16, color=MUTED)
notes(s, "This matters for comparing with headline depth: the sonar depth depends on the tilt, not on where the fish are.")

# ---------------------------------------------------------------------------------------------- 13 ST90 vs SU90
s = prs.slides.add_slide(BLANK)
title(s, "ST90 and SU90 side by side")
cd = CategoryChartData()
cd.categories = ["ST90", "SU90"]
cd.add_series("Usable schools", (4034, 2472))
gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(0.6), Inches(1.5), Inches(4.2), Inches(4.6), cd)
ch = gf.chart
ch.has_legend = False
ch.has_title = True
ch.chart_title.text_frame.text = "Usable schools, same hours"
ch.plots[0].has_data_labels = True
ch.plots[0].data_labels.number_format, ch.plots[0].data_labels.number_format_is_linked = "#,##0", False
ch.plots[0].data_labels.position = XL_LABEL_POSITION.OUTSIDE_END
ch.plots[0].gap_width = 60
ch.value_axis.visible = False
ch.value_axis.has_major_gridlines = False
cd2 = CategoryChartData()
cd2.categories = ["ST90", "SU90"]
cd2.add_series("Tracks removed", (0.11, 0.29))
gf2 = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(5.0), Inches(1.5), Inches(3.6), Inches(4.6), cd2)
ch2 = gf2.chart
ch2.has_legend = False
ch2.has_title = True
ch2.chart_title.text_frame.text = "Tracks removed as unrealistic"
ch2.plots[0].has_data_labels = True
ch2.plots[0].data_labels.number_format, ch2.plots[0].data_labels.number_format_is_linked = "0%", False
ch2.plots[0].data_labels.position = XL_LABEL_POSITION.OUTSIDE_END
ch2.plots[0].gap_width = 60
ch2.value_axis.visible = False
ch2.value_axis.has_major_gridlines = False
ch2.value_axis.maximum_scale = 0.35
for c in (ch, ch2):
    c.chart_title.text_frame.paragraphs[0].runs[0].font.size = Pt(16)
    c.chart_title.text_frame.paragraphs[0].runs[0].font.bold = True
    c.chart_title.text_frame.paragraphs[0].runs[0].font.color.rgb = NAVY
    c.category_axis.tick_labels.font.size = Pt(16)
    c.plots[0].data_labels.font.size = Pt(16)
    c.plots[0].data_labels.font.bold = True
    for i, pt in enumerate(c.plots[0].series[0].points):
        pt.format.fill.solid()
        pt.format.fill.fore_color.rgb = ORANGE if i == 0 else TEAL
text(s, 9.0, 1.6, 3.9, 5.0, [
    [("ST90: ", {"bold": True, "color": ORANGE}), ("more usable schools, followed for longer and a bit more precisely.", {})],
    [("Both agree: ", {"bold": True, "color": NAVY}), ("the same west-south-west movement of the herring.", {})],
    [("Not a fully fair test: ", {"bold": True, "color": NAVY}), ("the ST90 was tilted at 4°, the SU90 at 9°, so they looked at different water.", {})]],
     size=17, space_after=14)
caption(s, "June–July 2026, hours when both sonars recorded. Only 396 schools were followed by both sonars at once, mostly because of the different tilts.")
notes(s, "Good news either way: the SU90, which is common on other vessels, tells the same story as the ST90.")

# ---------------------------------------------------------------------------------------------- 14 next
s = prs.slides.add_slide(BLANK)
title(s, "How you could help next season")
asks = [("A day with matched tilts", "Both sonars at the same tilt for a day or two, so we can calibrate one against the other"),
        ("Clocks on UTC", "The ST90 clock ran about one hour ahead of the SU90 and the catch log"),
        ("Keep recording", "Autumn and other seasons, to see if the westward movement comes back each year"),
        ("Your observations", "What you saw on the sonar and on deck, to check our results against")]
for i, (h, b) in enumerate(asks):
    x = 0.6 + (i % 2) * 6.2
    y = 1.6 + (i // 2) * 2.55
    shape(s, MSO_SHAPE.ROUNDED_RECTANGLE, x, y, 5.9, 2.25, TINT).adjustments[0] = 0.08
    badge(s, x + 0.3, y + 0.3, str(i + 1))
    text(s, x + 1.1, y + 0.3, 4.5, 0.5, h, size=20, bold=True, color=NAVY)
    text(s, x + 1.1, y + 0.85, 4.5, 1.3, b, size=16, color=MUTED)
notes(s, "These are requests, not obligations. The tilt experiment is the most valuable and only needs a day or two.")

# ---------------------------------------------------------------------------------------------- 15 close
s = prs.slides.add_slide(BLANK)
dark_bg(s)
text(s, 0.8, 0.8, 11, 1.0, "In short", size=40, font=HEAD, bold=True, color=WHITE)
take = [("The sonar can follow fish schools", "after cleaning out the ship's own echoes and smoothing the tracks."),
        ("Herring moved about 5 nm a day west-south-west", "by day, while the current was weak and ran south-east: they were swimming, possibly towards the spawning grounds."),
        ("Mackerel in 2025 were mostly carried by the current.", ""),
        ("School counts did not predict catch;", "the sonar tells us more about how fish behave than about how much you will catch.")]
for i, (h, b) in enumerate(take):
    y = 2.0 + i * 1.1
    badge(s, 0.8, y, str(i + 1))
    text(s, 1.6, y + 0.03, 10.8, 0.9, [[(h + " ", {"bold": True, "color": WHITE}), (b, {"color": RGBColor(0xCF, 0xDD, 0xE8)})]], size=20)
text(s, 0.8, 6.5, 11, 0.6, "Thank you. Questions?", size=24, bold=True, color=ORANGE)
notes(s, "Close with thanks and invite questions; ask what he recognises from his own experience.")

os.makedirs(OUT, exist_ok=True)
f = os.path.join(OUT, "SCH302_sonar_results_for_skipper.pptx")
prs.save(f)
print("saved", f)
