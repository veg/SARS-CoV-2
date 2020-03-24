/*
add download button on top left
*/

document.body.style.border = "3px dashed blue";

let dl_id = 'scraper_button';
let dl_button = document.getElementById(dl_id);

let sequence_json = {};
let records = [];
let current_record = 0;

var observer = new MutationObserver(function(mutations) {
  mutations.forEach(function(mutation) {
    if (!mutation.addedNodes) return;

    for (var i = 0; i < mutation.addedNodes.length; i++) {
      var node = mutation.addedNodes[i];      
      if (node.nodeName == "IFRAME") {
        node.onload = override_onload;  
      }
    }
    
    for (var i = 0; i < mutation.removedNodes.length; i++) {
      var node = mutation.removedNodes[i];      
      if (node.nodeName == "IFRAME") {
        current_record ++;
        if (current_record >= records.length) {
            dl_button.innerHTML = "Download " + Object.keys(sequence_json).length + " records to file";
        } else {
            records[current_record].click();
        }
      }
    }
    
  })
});

function downloadObjectAsJson(exportObj, exportName){
    var dataStr = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(exportObj));
    var downloadAnchorNode = document.createElement('a');
    downloadAnchorNode.setAttribute("href",     dataStr);
    downloadAnchorNode.setAttribute("download", exportName + ".json");
    document.body.appendChild(downloadAnchorNode); // required for firefox
    downloadAnchorNode.click();
    downloadAnchorNode.remove();
}

function onStartedDownload(id) {
  console.log(`Started downloading: ${id}`);
}

function onFailed(error) {
  console.log(`Download failed: ${error}`);
}

function scrape_records (event) {
    if (records.length == 0) {
        records = document.body.querySelectorAll("tr.yui-dt-rec");
        if (records.length) {
            current_record = 0;
            records[0].click();
            observer.observe(document.body, {
                childList: true
              , subtree: false
              , attributes: false
              , characterData: true
            });
        }
    } else {
        if (current_record >= records.length) {
            downloadObjectAsJson (sequence_json, 'gisaid');
            records = [];
            dl_button.innerHTML =  "Scrape records from this page";
        } else {
            observer.disconnect();
            current_record = records.length;
            dl_button.innerHTML = "Scrape records from this page";
        }
    }
}

function override_onload (event) {
    //console.log ("IFRAME loaded", event);
    dl_button.innerHTML = "Downloading sequence " + (current_record + 1);
    let fields = event.target.contentDocument.querySelectorAll("tr");
    let sequence_record = {};
    for (let i = 0; i < fields.length; i++) {
        let cells = fields[i].querySelectorAll("td");
        if (cells.length == 2) {
            let value = cells[1].innerText;
            if (value.length) {
                let label = cells[0].innerText.split (':')[0];
                sequence_record[label] = value;
            }
            
        }
    }
    if ("Accession ID" in sequence_record) {
        let seq_data = event.target.contentDocument.querySelector ("pre");
        if (seq_data) {
            sequence_record ["FASTA"] = seq_data.textContent;
            sequence_json [sequence_record["Accession ID"]] = sequence_record;
        }
    }
    
        
    
    dl_button.innerHTML = "Closing window";
    let closeme = event.target.contentDocument.querySelectorAll("button");
    if (closeme.length >= 2) {
        closeme[1].click();
    }
    dl_button.innerHTML = "Extracted " + Object.keys(sequence_json).length + " records";
}


if (!dl_button) {
    dl_button = document.body.appendChild( document.createElement('button'));
    dl_button.id = dl_id;
}

dl_button.innerHTML = "Scrape records from this page";
dl_button.onclick = scrape_records;







