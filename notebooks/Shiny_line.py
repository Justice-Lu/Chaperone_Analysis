import os 
import numpy as np 
import pandas as pd 
import copy
import matplotlib.pyplot as pl
# following line ensures the pl plots inline 
import plotly.express as px
import plotly.graph_objects as go

from importlib import reload
import GE_functions

ge_normalized = pd.read_csv('../expression_csv/ge_normalized_GSE151346_MOE_ALL_OlfrSum.csv', 
                            index_col=0)



from shiny import App, reactive, render, ui

custom_label = "Some Counter modified "

counter_ui = ui.div(
    {"style": "border: 1px solid #ccc; border-radius: 5px; margin: 5px 0;"},
    ui.h2("This is " + custom_label),
    ui.input_action_button(id="button", label=custom_label),
    ui.output_text_verbatim(id="out"),
)

def counter_server(input, output, session, starting_value = 0):
    count = reactive.Value(starting_value)

    @reactive.Effect
    @reactive.event(input.button)
    def _():
        count.set(count() + 1)

    @output
    @render.text
    def out():
        return f"Click count is {count()}"

app = App(counter_ui, counter_server)