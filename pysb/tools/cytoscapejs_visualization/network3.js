/**
 * Created by dionisio on 2/22/17.
 */
$(function(){ // on dom ready

    var path = '/visualization_cytoscapejs/visualization_batch/graph_data2.js';
    var json_data = '';

    $.ajax({
        url: path,
        dataType: 'json',
        success: function(data){
            json_data = data;
        },
        async: false // turn into sync, so it can deal with global variable
    });

    var tspan = json_data.data.tspan;
    var model_name = json_data.data.name

    document.getElementById('myHeader').innerHTML = model_name;
    document.getElementById('range').max = tspan.length;

    $('#cy').cytoscape({
        boxSelectionEnabled: false,
        autounselectify: true,

        layout: {
            name: 'preset'
        },

        style: cytoscape.stylesheet()
            .selector('node')
            .css({
                'label': 'data(label)',
                'pie-size': '80%',
                'pie-1-background-color': '#916712',
                'pie-1-background-size': '100',
                'pie-2-background-color': '#dddcd4',
                'pie-2-background-size': '100',
            })

            .selector('edge')
            .css({
                'curve-style': 'bezier',
                // 'width': 'data(width)',
                'target-arrow-shape': 'data(target_arrow_shape)',
                'source-arrow-shape': 'data(source_arrow_shape)',
            }),

        elements: json_data.elements,

        ready: function(){
            window.cy = this;

            // giddy up
        }



    });


        // console.log(`data(edge_size_t${t})`)
    // for (let t = 0; t < 99; t++) {
    //     cy.batch(function () {
    //         cy.edges().forEach(function (n) {
    //             var c = n.data(`edge_color_t${t}`);
    //             n.style({'line-color': c})
    //
    //         })
    //     });
    //     console.log(cy.edges()[0].style())
    //
    // }

    var idx = 0;
    var t = 0;

    cy.elements().qtip({
        content: ' ',
        position: {
            my: 'top center',
            at: 'bottom center'
        },
        style: {
            classes: 'qtip-bootstrap',
            tip: {
                width: 8,
                height: 4
            }
        }
    });

    function update_graph(t){
        cy.batch(function () {
            cy.edges().forEach(function (e) {
                var c = e.data(`edge_color_t${t}`);
                var s = e.data(`edge_size_t${t}`);
                var a = e.data(`edge_qtip_t${t}`);
                e.animate({
                    style: {'line-color': c,
                        'target-arrow-color': c,
                        'source-arrow-color': c,
                        'width': s},
                    duration: 1000
                });
                e.qtip('api').set('content.text', a.toString());
                // n.animate({style: {'width': s}})

            });
            cy.nodes().forEach(function(n){
                var p = n.data(`rel_value_t${t}`);
                var q = n.data(`abs_value_t${t}`);
                n.animate({
                    style: {'pie-1-background-size': p},
                    duration: 1000
                });
                n.qtip('api').set('content.text', q.toString())
            })
        })
    }

    var animating = (function($){
        var GuiPause = $('#Pause');
        var GuiResume = $('#Resume').hide();

        var Running = true;
        var advancing = function(){
            if (idx < tspan.length){
                update_graph(idx)
                ;

                t = setTimeout(advancing, 1000)
            }
            if (Running){
                idx++;
            }
            // else{
            //     console.log('not running')
            // }
        };

        var Pause = function() {
            Running = false;
            GuiPause.hide();
            GuiResume.show();
        };

        var Resume = function() {
            Running = true;
            GuiPause.show();
            GuiResume.hide();
        };

        var Start = function(){
            advancing();
        };
        return {
            Pause: Pause,
            Resume: Resume,
            Start: Start
        };
    })((jQuery));

    // jQuery('#pause').on('click',animating.Pause);
    // jQuery('#resume').on('click',animating.Resume);

    var playButton = document.getElementById('Play');
    playButton.addEventListener('click', function(){
        animating.Start();
    });

    var pauseButton = document.getElementById('Pause');
    pauseButton.addEventListener('click', function() {
        animating.Pause();
    });

    var resumeButton = document.getElementById('Resume');
    resumeButton.addEventListener('click', function() {
        animating.Resume();
    });

    var rangeInput = document.getElementById("range");

    rangeInput.addEventListener('mouseup', function() {
        update_graph(this.value)
    });

    var range = document.getElementById('range'),
        text = document.getElementById('textInput');

    range.onchange = function() {
        text.value = tspan[this.value].toFixed(2);
        idx = this.value
    }

}); // on dom ready