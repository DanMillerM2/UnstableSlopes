# Quantiles
Dan
2025-03-26

Even after alignment of two sequential DTMs, there will be persistent
errors in the measured elevation differences. These errors arise
primarily from wide ground-return spacing and reflections from
vegetation interpreted as ground returns. This is shown in
<a href="#fig-errors" class="quarto-xref">Figure 1</a> for a small area
from the 2006-2017 overlap examined in the Coregistration document. The
left panel shows elevation differences with a linear stretch from -0.3m
(blue) to +0.3m (red). Yellow indicates a difference near zero. Over the
low-gradient, low-relief, unforested fields and pastures, elevation
differences are near zero with little variability. On the forested
hillslopes, the speckled appearance indicates much greater variability
in elevations between the two DTMs.

<div id="fig-errors">

<img src="images/clipboard-113032948.png" data-fig-align="center" />

Figure 1: Noise in the elevation differences increases in areas with
tree canopy.

</div>

We want to differentiate ground-surface elevation changes from noise in
the two DTMs. One approach is to determine the minimum difference above
which differences between the two DTMs can reliably be interpreted as
changes in elevation. This minimum difference is referred to as the
level of detection (LoD). Error statistics compiled for each DTM based
on comparison of ground control points with lidar ground-return
elevations can be used to estimate a theoretical LoD based on concepts
of error propagation (Brasington, Langham, and Rumsby 2003; Lane,
Westaway, and Murray Hicks 2003). Using the standard deviation of
elevation errors estimated for each DTM ($\sigma_{2006}$ for the 2006
lidar and $\sigma_{2017}$ for the 2007 lidar), the standard deviation of
differences in elevations between the DTMs ($\sigma_{\Delta z}$) is
estimated as

$$
\sigma_{\Delta z} = \sqrt{\sigma_{2006}^2 + \sigma_{2017}^2}
$$

For the 2006 lidar, $\sigma_{2006}$ is reported as 0.06m. For the 2017
lidar, $\sigma_{2017}$ is reported as 0.04 m in unvegetated areas and
0.09m in vegetated areas. Using the vegetated value, this gives an
estimated $\sigma_{\Delta z}$ of about 0.10 m. Assuming that errors
follow a normal distribution centered about zero, we might then assume
with 95% confidence, using two $\sigma_{\Delta Z}$, that differences
greater than 0.2m indicate actual changes in elevation. The speckled
pattern over forested areas visible in the left panel of
<a href="#fig-errors" class="quarto-xref">Figure 1</a> above, however,
indicates apparently random variations in $\Delta z$ considerably larger
than 0.2m over distances of ten meters or less. In those areas, an LoD
of 0.2m would those red and blue speckles. In contrast, over the
unforested fields and pastures, variations in $\Delta z$ are
considerably less than 0.2m and an LoD less than 0.2m may be
appropriate.

What is the source of this variability?

![](images/clipboard-1515199096.png)

We might also use the results reported in the Coregistration document,
where we used the interquartile range in $\Delta z$ measured over the
entire study area to estimate the limits of detectability. We found a
slope-dependent relationship, from a limit of about 0.1m on flat terrain
up to a limit of about 0.25m on steep terrain. Given, however, the
increase in uncertainty in elevation differences under forest canopy,
this may simply reflect a trend for fields and pastures to occur in
flatter areas and forests on hillslopes. We need a method to measure
local variability in $\Delta z$ with which to estimate our level of
certainty that a difference in elevation between two DTMs represents an
actual change in elevation. We can expand on the use of quartiles in
measured $\Delta z$ values that we developed to identify stable terrain
for coregistration.

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-brasington2003" class="csl-entry">

Brasington, James, Joe Langham, and Barbara Rumsby. 2003.
“Methodological Sensitivity of Morphometric Estimates of Coarse Fluvial
Sediment Transport.” *Geomorphology* 53 (3-4): 299–316.
<https://doi.org/10.1016/S0169-555X(02)00320-3>.

</div>

<div id="ref-lane2003" class="csl-entry">

Lane, Stuart N., Richard M. Westaway, and D. Murray Hicks. 2003.
“Estimation of Erosion and Deposition Volumes in a Large, Gravel-Bed,
Braided River Using Synoptic Remote Sensing.” *Earth Surface Processes
and Landforms: The Journal of the British Geomorphological Research
Group* 28 (3): 249–71.

</div>

</div>
