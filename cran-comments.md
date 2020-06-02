## Test environments
* local R installation, R 4.0.0
* win-builder (devel)

## R CMD check results
There were no ERRORs or WARNINGs

### local R installation (Mac), R 4.0.0
0 errors | 0 warnings | 1 note

* This is a new release.

### win-builder (devel)
0 errors | 0 warnings | 1 note

R Under development (unstable) (2020-05-22 r78545)

## Reply

### \dontrun{} in plot_omega() function
The code in \dontrun{} cannot be executed.
The first one (nawt()) can be executed but the second one (plot_omega()) 
cannot because it needs the value of "omega" in the "object" argument but
nawt(method = "both", twostep = FALSE) does not return omega.


