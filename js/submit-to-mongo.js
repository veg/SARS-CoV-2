const _ = require("lodash"),
  readline = require('readline'),
  fs = require("fs"),
  util = require("util"),
  commander = require("commander"),
  moment = require("moment"),
  mongodb = require("mongodb");

// example gisaid record
//{
//  "b": "EPI_ISL_402119",
//  "c": "hCoV-19/Wuhan/IVDC-HB-01/2019",
//  "d": "Virus Isolate, Passage 1",
//  "e": "EPI_ISL_402119",
//  "f": "2019-12-30",
//  "g": "2020-01-10",
//  "h": null,
//  "i": 29891,
//  "j": "Human",
//  "k": "Asia / China / Hubei / Wuhan",
//  "l": "National Institute for Viral Disease Control and Prevention, China CDC",
//  "m": "National Institute for Viral Disease Control and Prevention, China CDC"
//},

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
let metadata = JSON.parse(fs.readFileSync(filename));
let records = metadata.records;

// remove duplicate seqs
console.log("Got " + records.length + " metadata records");

let adaptedRecords = _.map(records, record => {

  // switch to just keep what was there.
  let collectionDate = record.g;
  let submissionDate = record.h;
  let locs = _.split(record.l, "/");

  let loc = {
    subregion: _.trim(locs[0]),
    country: _.trim(locs[1]),
    state: _.trim(locs[2]),
    locality: _.trim(locs[3])
  };

  let acc = _.toLower(record.b);

  let adaptedRecord = {
    address: null,
    age: null,
    assembly: null,
    authors: null,
    collected: collectionDate,
    coverage: null,
    gender: null,
    host: record.k,
    id: acc,
    lab: record.m,
    location: loc,
    name: record.d,
    passage: record.e,
    seqLength: record.j,
    submitted: submissionDate,
    submitter: record.n,
    technology: null,
    type: null
  };

  return adaptedRecord;

});

const MongoClient = mongodb.MongoClient;
const url = 'mongodb://129.32.209.134:27017';
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

	collection.insertMany(adaptedRecords, {'ordered' : false }, (err, result) => {
		if (err) {
			console.log(err);
		}

  	client.close();
    process.exit(1)

	});


});

