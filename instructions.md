# Instructions on how to run ReTime to predict Retention Times for Lipids
## Basics
For this to work you will need to upload a datafile of lipids and their retention times.
Input should have "Train" and "Test" labelled in a csv as below

-file is in csv format
-file has a column "Lipid ID" (without the quotations) which specifies all the lipid names in conventional format (e.g. Cer (d18:1/20:1)).
-file has a column "Retention Time" (without the quotations) which specifies the user recorded retention times
-file has a column "Type" (without the quotations) with values either "Train" or "Test" (without the quotations) - note any -record with "Train" must have a corresponding retention time


Run your model and view the predicted Retention times for the Lipids you are looking to predict


###### University of Sydney