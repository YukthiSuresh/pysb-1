import pydot
import pysb.viz.graphviz as viz
from pysb.examples.tyson_oscillator import model


def test_render_reaction():
    rxn_graph = viz.render_reactions(model)
    assert isinstance(rxn_graph, pydot.Dot)
    rxn_graph = viz.run_render_reactions(model, to_string=True)
    assert isinstance(rxn_graph, (str, bytes))


def test_render_species():
    species_graph = viz.run_render_species(model)
    assert isinstance(species_graph, pydot.Dot)
    species_graph = viz.run_render_species(model, None, to_string=True)
    assert isinstance(species_graph, (str, bytes))


if __name__ == '__main__':
    test_render_species()
    test_render_reaction()
