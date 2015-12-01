 encode site_side, generate(ss)
 encode eventual_mpa, generate(will_be_mpa)

generate arg= fishable* mpa

xi: tobit logdensity fishable mpa arg ss mpa_area, ll(-2.3025851)

predict predicted2

plot logdensity predicted2

