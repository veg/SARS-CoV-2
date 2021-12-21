import puppeteer from "https://deno.land/x/puppeteer@9.0.0/mod.ts";
import { cheerio } from "https://deno.land/x/cheerio@1.0.4/mod.ts";
import * as path from "https://deno.land/std@0.61.0/path/mod.ts";
import * as R from "https://deno.land/x/ramda@v0.27.2/mod.ts";
import  __ from 'https://deno.land/x/dirname/mod.ts';
const { __filename, __dirname } = __(import.meta);

const url = 'https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/index.html';

/**
 * write.ts
 */
function writeJson(path: string, data: object): string {
  try {
    Deno.writeTextFileSync(path, JSON.stringify(data));

    return "Written to " + path;
  } catch (e) {
    return e.message;
  }
}

try {
    const browser = await puppeteer.launch({
      headless: true,
      args: [
        '--no-sandbox',
        '--disable-setuid-sandbox',
      ]
    });
    const page = await browser.newPage();
    await page.goto(url);

    const html = await page.content();

    const $ = cheerio.load(html);


    let cladeNames: any[] = []

    $('table').each((index, element) => {
      $(element).find('tr').each((rowIndex, rowElement) => {
        if(index < 2) {
          console.log($(rowElement).children('td').eq(1).text());
          let text = $(rowElement).children('td').eq(1).text().replace('#', '').trim();
          cladeNames.push(text);
        } else {
          let text = $(rowElement).children('td').eq(0).text().trim();
          console.log(text);
          cladeNames.push(text);
        }
      })
      return cladeNames
    });

    let noEmpty = R.reject(R.isEmpty, cladeNames);
    let removeCA = R.reject((d: string) => R.includes("B.1.427", d), noEmpty);
    let removeMal = R.reject((d: string) => R.includes("_", d) || R.includes("#", d) || R.match(/[0-9]/, R.last(d)).length == 0, removeCA);

    let clades = {
      'clades': removeMal
    }

    writeJson(__dirname + '/../airflow/libs/voc.json', clades)
    Deno.exit(0);


} catch(error) {
    console.log(error);
}
