#!/usr/bin/env python
"""
Example
-------

    $ python -m pysb.viz.graphviz render_reaction mymodel.py > mymodel_reactions.dot
    $ python -m pysb.viz.graphviz  render_species mymodel.py > mymodel_species.dot
    $ python -m pysb.viz.graphviz render_reaction mymodel.py model_reactions.png
    $ python -m pysb.viz.graphviz  render_species mymodel.py model_species.png

Renders the reactions produced by a model into the "dot" graph format which can
be visualized with Graphviz.

To create a PDF from the .dot file, use the "dot" command from Graphviz::

    dot mymodel_reactions.dot -T pdf -O

This will create mymodel.dot.pdf. You can also change the "dot" command to one
of the other Graphviz drawing tools for a different type of layout. Note that
you can pipe the output of render_reactions straight into Graphviz without
creating an intermediate .dot file, which is especially helpful if you are
making continuous changes to the model and need to visualize your changes
repeatedly::

    python -m pysb.viz.graphviz render_reaction mymodel.py | dot -T pdf -o mymodel_reactions.pdf

Note that some PDF viewers will auto-reload a changed PDF, so you may not even
need to manually reopen it every time you rerun the tool.

Output for Robertson example model
==================================

The Robertson example model (in ``pysb.examples.robertson``) contains
the following three reactions:

* A -> B
* B + B -> B + C
* C + B -> C + A

The reaction network diagram for this system as generated by this module and
rendered using ``dot`` is shown below:

.. image:: robertson_reactions.png
   :align: center
   :alt: Reaction network for pysb.examples.robertson

Circular nodes (``r0``, ``r1`` and ``r2``) indicate reactions; square nodes
(``A()``, ``B()`` and ``C()``) indicate species. Incoming arrows from a species
node to a reaction node indicate that the species is a reactant; outgoing
arrows from a reaction node to a species node indicate that the species is a
product. A hollow diamond-tipped arrow from a species to a reaction indicates
that the species is involved as both a reactant and a product, i.e., it serves
as a "modifier" (enzyme or catalyst).

"""

from __future__ import print_function
import sys
import os
import re
import pysb.bng
import pydot


# Alias basestring under Python 3 for forwards compatibility
try:
    basestring
except NameError:
    basestring = str


def run_render_species(model, save_name=None):
    """
    Render the species from a model into the "dot" graph format.

    Parameters
    ----------
    model : pysb.core.Model
        The model to render.
    save_name : str, optional
        Will render and save an image to png with name provided
    Returns
    -------
    string
        The dot format output.
    """

    pysb.bng.generate_equations(model)
    # return render_species_as_dot(model.species, model.name)
    return render_species(model.species, model.name, save_name
                          ).create(prog='dot', format='dot')


def run_render_reactions(model, save_name=None):
    pysb.bng.generate_equations(model)
    # return render_species_as_dot(model.species, model.name)
    return render_reactions(model, save_name).create(prog='dot', format='dot')


def render_species(species_list, graph_name="", save_name=None):
    graph = pydot.Dot(graph_name="{} species".format(graph_name),
                      graph_type='digraph', rankdir="LR", fontname='Arial',
                      dpi=200)

    graph.set_edge_defaults(fontname='Arial', fontsize=8)
    for si, cp in enumerate(species_list):
        sgraph_name = 'cluster_s%d' % si
        cp_label = re.sub(r'% ', '%<br align="left"/>',
                          str(cp)) + '<br align="left"/>'
        sgraph_label = '<<font point-size="10" color="blue">s%d</font>' \
                       '<br align="left"/>' \
                       '<font face="Consolas" point-size="6">%s</font>>' % \
                       (si, cp_label)
        subgraph = pydot.Cluster(graph_name=sgraph_name, label=sgraph_label,
                                 color="gray75", sortv=sgraph_name)

        bonds = {}
        for mi, mp in enumerate(cp.monomer_patterns):
            monomer_node = '%s_%d' % (sgraph_name, mi)
            monomer_label = '<<table border="0" cellborder="1" cellspacing="0">'
            monomer_label += '<tr><td bgcolor="#a0ffa0"><b>%s</b></td></tr>' % \
                             mp.monomer.name
            for site in mp.monomer.sites:
                site_state = None
                cond = mp.site_conditions[site]
                if isinstance(cond, basestring):
                    site_state = cond
                elif isinstance(cond, tuple):
                    site_state = cond[0]
                site_label = site
                if site_state is not None:
                    site_label += '=<font color="purple">%s</font>' % site_state
                monomer_label += '<tr><td port="%s">%s</td></tr>' % \
                                 (site, site_label)
            monomer_label += '</table>>'
            for site, value in mp.site_conditions.items():
                site_bonds = []  # list of bond numbers
                if isinstance(value, int):
                    site_bonds.append(value)
                elif isinstance(value, tuple):
                    site_bonds.append(value[1])
                elif isinstance(value, list):
                    site_bonds += value
                for b in site_bonds:
                    bonds.setdefault(b, []).append((monomer_node, site))

            node = pydot.Node(monomer_node, label=monomer_label, shape="none",
                              fontname="Arial", fontsize=8)
            subgraph.add_node(node)
        for bi, sites in bonds.items():
            node_names, port_names = list(zip(*sites))
            edge = pydot.Edge(node_names[0], node_names[1],
                              ltail=port_names[0], lhead=port_names[1],
                              label=str(bi)
                              )
            subgraph.add_edge(edge)
        graph.add_subgraph(subgraph)
    if save_name is not None:
        _save(graph, save_name)
    return graph


def render_reactions(model, save_name=None):
    """
    Render the reactions produced by a model into the "dot" graph format.

    Parameters
    ----------
    model : pysb.core.Model
        The model to render.
    save_name : str, optional
        Will render and save an image to png with name provided
    Returns
    -------
    string
        The dot format output.
    """

    pysb.bng.generate_equations(model)

    graph = pydot.Dot(directed=True, rankdir="LR", dpi=300)
    ic_species = [cp for cp, parameter in model.initial_conditions]
    for i, cp in enumerate(model.species):
        species_node = 's%d' % i
        slabel = re.sub(r'% ', r'%\\l', str(cp))
        slabel += '\\l'
        color = "#ccffcc"
        # color species with an initial condition differently
        if len([s for s in ic_species if s.is_equivalent_to(cp)]):
            color = "#aaffff"

        graph.add_node(
            pydot.Node(species_node,
                       label=slabel,
                       shape="Mrecord",
                       fillcolor=color, style="filled", color="transparent",
                       fontsize="12",
                       margin="0.06,0")
        )
    for i, reaction in enumerate(model.reactions_bidirectional):
        reaction_node = 'r%d' % i
        graph.add_node(
            pydot.Node(reaction_node,
                       label=reaction_node,
                       shape="circle",
                       fillcolor="lightgray", style="filled", color="transparent",
                       fontsize="12",
                       width=".3", height=".3", margin="0.06,0"))
        reactants = set(reaction['reactants'])
        products = set(reaction['products'])
        modifiers = reactants & products
        reactants = reactants - modifiers
        products = products - modifiers
        attr_reversible = {'dir': 'both', 'arrowtail': 'empty'} if reaction['reversible'] else {}
        for s in reactants:
            r_link(graph, s, i, **attr_reversible)
        for s in products:
            r_link(graph, s, i, _flip=True, **attr_reversible)
        for s in modifiers:
            r_link(graph, s, i, arrowhead="odiamond")
    if save_name is not None:
        _save(graph, save_name)
    return graph


def _save(graph, save_name):
    if len(save_name.split('.')) == 2:
        prefix, frmt = save_name.split('.')
        if frmt in graph.formats:
            graph.write(save_name, prog='dot', format=frmt)
    else:
        print("Unknown format, assuming pdf")
        graph.write_pdf('{}.pdf'.format(save_name))


def r_link(graph, s, r, **attrs):
    nodes = ('s%d' % s, 'r%d' % r)
    if attrs.get('_flip'):
        del attrs['_flip']
        nodes = list(reversed(nodes))
    attrs.setdefault('arrowhead', 'normal')
    graph.add_edge(pydot.Edge(nodes[0], nodes[1], **attrs))


usage = __doc__
usage = usage[1:]  # strip leading newline

if __name__ == '__main__':
    # sanity checks on filename
    if len(sys.argv) <= 2:
        print(usage, end=' ')
        exit()
    model_filename = sys.argv[2]
    if len(sys.argv) == 4:
        save_name = sys.argv[3]
    else:
        save_name = None
    print(save_name)
    print(sys.argv)
    if not os.path.exists(model_filename):
        raise Exception("File '%s' doesn't exist" % model_filename)
    if not re.search(r'\.py$', model_filename):
        raise Exception("File '%s' is not a .py file" % model_filename)
    sys.path.insert(0, os.path.dirname(model_filename))
    model_name = re.sub(r'\.py$', '', os.path.basename(model_filename))
    # import it
    try:
        # FIXME if the model has the same name as some other "real" module
        # which we use, there will be trouble
        # (use the imp package and import as some safe name?)
        model_module = __import__(model_name)
    except Exception as e:
        print("Error in model script:\n")
        raise
    # grab the 'model' variable from the module
    try:
        model = model_module.__dict__['model']
    except KeyError:
        raise Exception("File '%s' isn't a model file" % model_filename)

    if sys.argv[1] == 'render_reactions':
        print(run_render_reactions(model, save_name=save_name))
    elif sys.argv[1] == 'render_species':
        print(run_render_species(model, save_name=save_name))
    else:
        print(usage, end=' ')
        exit()
