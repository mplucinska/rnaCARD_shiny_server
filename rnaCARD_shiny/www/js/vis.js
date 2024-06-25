
$(document).ready(function() {
  
  // create a click handler which listens for a click on the element with id equal to RStudio
  $("#submit").on("click", function(){
  
    Shiny.addCustomMessageHandler("mymessage", function(message) {
        var container = new fornac.FornaContainer("#rna1", { 'friction': 0 ,'layout':'naview','applyForce': true, 'allowPanningAndZooming': true});
    
        var options = {'structure': message.str1,
                       'sequence':  message.sequence,
                       'avoidOthers': true,
        };
        
        container.addRNA(options.structure, options);
        container.addCustomColorsText(message.color);
        container.setSize();
    });

    Shiny.addCustomMessageHandler("mymessage2", function(message) {
        var container2 = new fornac.FornaContainer("#rna2", { 'friction': 0 ,'layout':'naview','applyForce': true, 'allowPanningAndZooming': true});
    
       var options2 = {'structure': message.str2,
                       'sequence':  message.sequence,
                     'avoidOthers': true};
        container2.addRNA(options2.structure, options2);
        container2.addCustomColorsText(message.color);
        container2.setSize();
    });
    
    Shiny.addCustomMessageHandler("message_motif", function(message) {
        var container = new fornac.FornaContainer("#rna_m", { 'friction': 0 ,'layout':'naview','applyForce': true, 'allowPanningAndZooming': true});
        var options = {'structure': message.str_m,
                       'sequence':  message.seq_m,
                       'avoidOthers': true,
        };
        container.addRNA(options.structure, options);
        container.addCustomColorsText(message.col);
        container.setSize();
    });
    
        Shiny.addCustomMessageHandler("message_motif2", function(message) {
        var container = new fornac.FornaContainer("#rna_m2", { 'friction': 0 ,'layout':'naview','applyForce': true, 'allowPanningAndZooming': true});
        var options = {'structure': message.str_m,
                       'sequence':  message.seq_m,
                       'avoidOthers': true,
        };
        container.addRNA(options.structure, options);
        container.addCustomColorsText(message.col);
        container.setSize();
    });

    
        Shiny.addCustomMessageHandler("mymessage3", function(message) {
        var container = new fornac.FornaContainer("#rna3", { 'friction': 0 ,'layout':'naview','applyForce': true, 'allowPanningAndZooming': true});
    
        var options = {'structure': message.str1,
                       'sequence':  message.sequence,
                       'avoidOthers': true,
        };
        
        container.addRNA(options.structure, options);
        container.addCustomColorsText(message.color);
        container.setSize();
    });

    Shiny.addCustomMessageHandler("mymessage4", function(message) {
        var container = new fornac.FornaContainer("#rna4", { 'friction': 0 ,'layout':'naview','applyForce': true, 'allowPanningAndZooming': true});
    
       var options2 = {'structure': message.str2,
                       'sequence':  message.sequence,
                     'avoidOthers': true};
        container.addRNA(options2.structure, options2);
        container.addCustomColorsText(message.color);
        container.setSize();
    });
  });
});