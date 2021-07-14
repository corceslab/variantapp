from wtforms import (Form, TextField, IntegerField, SubmitField, validators)

class form_values(Form):
    """User entry form for entering cell type, chromosome, start, end, allele1, allele2"""
    cell_type = TextField("Enter a cell type:", validators=[validators.InputRequired()])
    chromosome = TextField("Enter the chromosome (chr#):", validators=[validators.InputRequired()])
    position = IntegerField("Enter the SNP position (0-based):", validators=[validators.InputRequired()])
    allele1 = TextField("Enter the effect base pair:", validators=[validators.InputRequired(), 
        validators.AnyOf(values=["A", "C", "G", "T"], message="allele1 must be A, C, G, or T")])
    allele2 = TextField("Enter the non-effect base pair:", validators=[validators.InputRequired(), 
        validators.AnyOf(values=["A", "C", "G", "T"], message="allele2 must be A, C, G, or T")])
    submit = SubmitField("Submit")

class form_rsID(Form):
    """User entry form for entering only an rsID"""
    cell_type = TextField("Enter a cell type:", validators=[validators.InputRequired()])
    rsID = TextField("Enter the variant's rsID (e.g. rs1237999):", validators=[validators.InputRequired()])
    submit = SubmitField("Submit")