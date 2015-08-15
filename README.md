# LODapp
[R], shinyapps, LoDapp - educational app for explaining the concept of limit of detection with a simulated chromatogram

### Description
This is a application that lets people working with LC-MS/MS systems to better understand the concept of limit of detection. The app creates a simulated chromatogram for which different parameters can be changed. The noise distibution and peak height distribution are also shown.

### Work to be done
At the moment the app in under development. I will be making the following changes to the app:
* Add the Critical Limit line to the noise distribution and show a calculation of the percent of peak height distribution that is below it (can be looked at as the p value of the test);
* Better descriptions will be added to the app;
* For some reason the updaing the mean value does not give a new chromatogram - this should be fixed;
* In the back end - the operations creating the chromatogram should be placed in one single function which is then used in the reactive({}).