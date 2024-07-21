### Prediction Model of Total Wealth
## Description
- I seek to find a model to predict the total wealth based on different feature variables about household information. The data I use comes from the 1991 Survey of Income and Program Participation (SIPP).

## Data
# Dependent Vairable:
- total wealth `tw`: the total market value of all the assets owned by a household. Total wealth equals net financial assets, including Individual Retirement Account (IRA) and 401(k) assets, plus housing equity plus the value of business, property, and motor vehicles”. 
# Feature variables 
- individual retirement account `ira`
- `e401`
- non-401k financial assets `nifa`
- income `inc`
- home mortgage `hmort`
- home value `hval`,
- home value minus home mortgage `hequity`
- educ: education `educ`
- `male`: 1 if male, 0 otherwise)
- `twoearn`: 1 if two earners in the household, 0 otherwise)
- education dummies: no highschool, high-school, some college, college,
- `age`
- family size `fsize`
- marital status `marr`: 1 if married, 0 otherwise
- There are a few terms that need to be clarified: individual retirement account (ira), e401. e401 plans are tax-deferred retirement savings accounts and it’s a company-sponsored account. They are offered by employers and only employees in companies offering such plans can participate while ira is an individual-created savings account. 
