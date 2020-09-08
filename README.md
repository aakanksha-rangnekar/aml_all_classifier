# aml_all_classifier
A classifier to distinguish between ALL and AML cancers. 

## Description 
Leukemia is a type of cancer that attack cells in the bone marrow that make blood. 
Acute lymphocytic lukemia(ALL) starts in the cells that become lymohocytes. Acute myelocytic lukemia begins in the early myeloid cells[1]. These two types of cancers have overlapping symptons making them harder to classify. The goal of the project is to develop a calssifier to distinguish betweeen AML and ALL using sequencing technlogy.
The data can be found under https://www.kaggle.com/crawford/gene-expression

## Development workflow 
The code starts by carrying out some exploratory data analysis, followed ny splitting the data into training and testing. Next the train data is normalized and selected features are used to build the model and to predict the train set. Once the model is finalized then the data is applied on the test set. A validation set can be used to then validate the model.


##### References
[1]https://www.webmd.com/cancer/lymphoma/leukemia-all-vs-aml#1

