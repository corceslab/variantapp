from types import ModuleType
from flask import Flask, request
from flask import render_template
from PIL import Image
import base64
import io
from markupsafe import Markup
import shutil

from utils.form import form_values, form_rsID
from utils.utils import generate_output_ranking

import tensorflow as tf
from tensorflow import compat


if __name__ == '__main__':
    cell_type = 'C24'
    rsID = 'rs691331, rs691342, rs691346'
    nc = '00'
    generate_output_ranking(cell_type, rsID, nc)
