const _ = require("lodash"),
  readline = require('readline'),
  d3 = require('d3-dsv'),
  fs = require("fs"),
  util = require("util"),
  commander = require("commander"),
  moment = require("moment"),
  mongodb = require("mongodb");

// example result record
//{
//  "epi_isl_402119" : {
//      "address": "National Institute for Viral Disease Control and Prevention, 155 Changbai Road, Changping District, Beijing 102206 China",
//      "age": "49",
//      "assembly": null,
//      "authors": "Wenjie Tan，Xiang Zhao，Wenling Wang，Xuejun Ma，Yongzhong Jiang，Roujian Lu, Ji Wang, Weimin Zhou，Peihua Niu，Peipei Liu，Faxian Zhan，Weifeng Shi，Baoying Huang，Jun Liu，Li Zhao，Yao Meng，Xiaozhou He，Fei Ye，Na Zhu，Yang Li，Jing Chen，Wenbo Xu，George F. Gao，Guizhen Wu",
//      "collected": "20191230",
//      "coverage": null,
//      "gender": "Female",
//      "host": "Human",
//      "id": "epi_isl_402119",
//      "lab": "National Institute for Viral Disease Control and Prevention, China CDC",
//      "location": {
//        "country": "China",
//        "locality": "Wuhan",
//        "state": "Hubei",
//        "subregion": "Asia"
//      },
//      "name": "hCoV-19/Wuhan/IVDC-HB-01/2019",
//      "passage": "Virus Isolate, Passage 1",
//      "submitted": "20200110",
//      "submitter": "Wenjie Tan",
//      "technology": null,
//      "type": "betacoronavirus"
//  }
//}

commander
  .arguments("<inputFile>", "Input filename")
  .action(cmd => {
    inputFile = cmd;
  })

commander
  .on("--help", function() {
    console.log("");
    console.log("Examples:");
    console.log("translate-to-master <file>");
  })
  .parse(process.argv);

let filename = inputFile;
//let metadata = d3.tsvParse(fs.readFileSync(filename));
let data = String(fs.readFileSync(filename));
let records = d3.tsvParse(data);

// remove duplicate seqs
console.log("Got " + records.length + " metadata records");


let adaptedRecords = _.map(records, record => {

  // switch to just keep what was there.
  let collectionDate = record['Collection date'];
  let submissionDate = record['Submission date'];

  let locationData = record['Location'].split('/');

  let loc = { 
    subregion : _.trim(locationData[0]), 
    supplemental : record['Additional location information'] 
  };

  if(locationData.length > 1) {
    loc['country'] = _.trim(locationData[1]);
  }

  if(locationData.length > 2) {
    loc['state'] = _.trim(locationData[2]);
  }

  if(locationData.length > 3) {
    loc['locality'] = _.trim(locationData[3]);
  }

  let acc = _.toLower(record['Accession ID']);
  let collectionDateType = '';

  try {
    collectionDateType = moment(collectionDate).toDate()
  } catch(e){}

  // Nextstrain_clade  pangolin_lineage  GISAID_clade
  let adaptedRecord = {
    address: null,
    genbank_accession: null,
    age: parseInt(record['Patient age']),
    assembly: null,
    authors: null,
    isComplete : record['Is complete?'],
    isGISAIDReference: record['Is reference?'],
    nContent : record['N-Content'],
    GCContent : record['GC-Content'],
    collected: collectionDateType,
    originalCollected: collectionDate,
    coverage: null,
    variant : record['Variant'],
    aaSubstitutions : record['AA Substitutions'],
    highCoverage: record['Is high coverage?'],
    lowCoverage: record['Is low coverage?'],
    length: parseInt(record['Sequence length']),
    gender: record['Gender'],
    sex: record['Gender'],
    host: record['Host'],
    id: acc,
    lab: null,
    originating_lab: null,
    location: loc,
    name: record['Accession ID'],
    passage: null,
    seqLength: parseInt(record['Sequence length']),
    submitted: submissionDate,
    submitted: moment(submissionDate).toDate(),
    originalSubmitted: submissionDate,
    submitter: null,
    submitting_lab: null,
    technology: null,
    type: record['Type'],
    pangolinLineage : record['Pango lineage'],
    pangolinVersion : record['Pangolin version'],
    gisaidClade : record['Clade']
  };

  return adaptedRecord;

});

const MongoClient = mongodb.MongoClient;
const url = 'mongodb://192.168.0.4:27017';
const dbName = 'gisaid';
const collectionName = 'records';

// Use the connect method to create a connection w/ the database
MongoClient.connect(url, (err, client) => {

  if (err) {
    throw err;
    return;
  }

  console.log('Database connection successful');

  // This objects holds the refrence to the db
  const db = client.db(dbName);
  const collection = db.collection(collectionName);

  collection.insertMany(adaptedRecords, {'ordered' : false, 'upsert': true }, (err, result) => {
    if (err) {
      console.log(err);
    }

    client.close();
    process.exit(0)

  });


});

