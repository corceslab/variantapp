from typing import Text
from wtforms import (Form, TextField, IntegerField, SelectField, SubmitField, validators, FieldList, FormField)

class form_entry(Form):
    pos = IntegerField("SNP Position")
    alleleA = TextField("Allele A:")
    alleleB = TextField("Allele B:")

class form_compound(Form):
    cell_type = SelectField(u'Cell Type', choices=[('C1 - Isocortical Excitatory'), ('C2 - Striatal Inhibitory (major)'), \
        ('C5 - Nigral Neurons (unclassified)'), ('C8 - OPCs'), \
        ('C13 - Astrocytes (unclassified)'), \
        ('C19 - Oligodendrocytes'), \
        ('C24 - Microglia'), \
        ('K562'), \
        ('GM12878')])
    chromosome = TextField("Chromosome (chr#):", validators=[validators.InputRequired()])
    position = IntegerField("Center Position:", validators=[validators.InputRequired()])
    num_used = IntegerField("Number of SNPs:", validators=[validators.InputRequired()])
    variants = FieldList(FormField(form_entry), min_entries = 5)
    submit = SubmitField("Submit")