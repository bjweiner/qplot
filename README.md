# qplot

qplot is an IDL routine for viewing spectra taken by the Hectospec instrument on
the MMT telescope and reduced with Richard Cool's HSRED pipeline. qplot allows
you to view the spectra in the output file made by HSRED, and displays the 
redshifts fitted by the HSRED reduce1d analysis. You can record redshift 
qualities and qplot will save the results to a text file for further use.

This version of qplot reads a catalog file that is made during the data reduction
using HSRED v 1.x. If you are using HSRED v2.0, you will probably have to alter the
location or format of the catalog file, which is just a text file with object IDs,
magnitude, etc for each fiber.

qplot was originally written by Casey Papovich, Christopher Willmer, and Ben Weiner,
and has been tuned up and released by Ben Weiner.

