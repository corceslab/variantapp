from typing import Text
from wtforms import (Form, TextField, IntegerField, SelectField, SubmitField, validators, FieldList, FormField)

class form_values(Form):
    """ User entry form for entering cell type, chromosome, start, end, allele1, allele2
    """
    
    cell_type = SelectField(u'Cell Type', choices=[('C1 - Isocortical Excitatory'), ('C2 - Striatal Inhibitory (major)'), \
        ('C5 - Nigral Neurons (unclassified)'), ('C8 - OPCs'), \
        ('C13 - Astrocytes (unclassified)'), \
        ('C19 - Oligodendrocytes'), \
        ('C24 - Microglia')])
    chromosome = TextField("Chromosome (chr#):", validators=[validators.InputRequired()])
    position = IntegerField("SNP position (0-based):", validators=[validators.InputRequired()])
    allele1 = TextField("Alternate Allele:", validators=[validators.InputRequired(), 
        validators.AnyOf(values=["A", "C", "G", "T"], message="alternate allele must be A, C, G, or T")])
    allele2 = TextField("Reference Allele:", validators=[validators.InputRequired(), 
        validators.AnyOf(values=["A", "C", "G", "T"], message="reference allele must be A, C, G, or T")])
    submit = SubmitField("Submit")

class form_rsID(Form):
    """User entry form for entering cell type (dropdown selection) and rsID (string query)
    """

    cell_type = SelectField(u'Cell Type', choices=[('C1 - Isocortical Excitatory'), ('C2 - Striatal Inhibitory (major)'), \
        ('C5 - Nigral Neurons (unclassified)'), ('C8 - OPCs'), \
        ('C13 - Astrocytes (unclassified)'), \
        ('C19 - Oligodendrocytes'), \
        ('C24 - Microglia')])
    rsID = TextField("Variant rsID List", validators=[validators.InputRequired()])
    submit = SubmitField("Submit")

class form_info(Form):
    rsIDs = TextField("Variant rsID List", validators=[validators.InputRequired()])
    submit = SubmitField("Submit")

class form_entry(Form):
    pos = IntegerField("SNP Position")
    alleleA = TextField("Allele A:")
    alleleB = TextField("Allele B:")

class form_compound(Form):
    cell_type = SelectField(u'Cell Type', choices=[('C1 - Isocortical Excitatory'), ('C2 - Striatal Inhibitory (major)'), \
        ('C5 - Nigral Neurons (unclassified)'), ('C8 - OPCs'), \
        ('C13 - Astrocytes (unclassified)'), \
        ('C19 - Oligodendrocytes'), \
        ('C24 - Microglia')])
    chromosome = TextField("Chromosome (chr#):", validators=[validators.InputRequired()])
    position = IntegerField("Center Position:", validators=[validators.InputRequired()])
    num_used = IntegerField("Number of Compound SNPs:", validators=[validators.InputRequired()])
    variants = FieldList(FormField(form_entry), min_entries = 10)
    submit = SubmitField("Submit")


# ('C1 - Isocortical Excitatory'), ('C2 - Striatal Inhibitory (major)'), ('C3 - Hippocampal Excitatory'), ('C4 - Hippocampal Excitatory'), \
#         ('C5 - Nigral Neurons (unclassified)'), ('C6 - Nigral Cells (unclassified)'), ('C7 - Neurons (unclassified)'), ('C8 - OPCs'), ('C9 - OPCs'), ('C10 - Nigral OPCs'), ('C11 - Isocortical Inhibitory'), \
#         ('C12 - Striatal Inhibitory (minor)'), ('C13 - Astrocytes (unclassified)'), ('C14 - Nigral Astrocytes'), ('C15 - Isocortical Astrocytes'), ('C16 - Striatal Astrocytes'), \
#         ('C17 - Astrocytes (unclassified)'), ('C18 - Potential Doublets'), ('C19 - Oligodendrocytes'), ('C20 - Oligodendrocytes'), ('C21 - Oligodendrocytes'), ('C22 - Oligodendrocytes'), \
#         ('C23 - Oligodendrocytes'), ('C24 - Microglia')