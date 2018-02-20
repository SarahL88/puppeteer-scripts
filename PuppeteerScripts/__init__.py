#!/usr/bin/env python
import datetime
import sys
import uuid
import requests
import json
from flask import request, Flask
app = Flask(__name__)

import modularcloning
import specificationview
import buildproject
import tecan_json_to_gwl

import tecan_to_autoprotocol

PUPPETEER_URL = 'http://128.31.25.45:8080/penguin-services-0.1.0/build'
CONSTELLATION_URL = 'http://localhost:8082/postSpecs'
CONCENTRATION_NG_UL = 25.0
CONCENTRATION_UNIT = 'NANOGRAMS_PER_MICROLITER'
VOLUME_UNIT = 'MICROLITERS'


@app.route("/get-tecan-instructions", methods=['POST'])
def get_tecan_instructions():
    authorid = str(uuid.uuid4())
    instanceid = str(uuid.uuid4())
    date = datetime.date.today()

    # Parse input files, populate 'repo' dict
    repo = modularcloning.make_repo(request.files, instanceid, authorid, date);

    # Get Constellation results, add to 'repo' dict
    specificationview.set_specification(repo, CONSTELLATION_URL, authorid, date);

    # Create output json
    puppeteer_input_json = buildproject.generate_build_request(repo, CONCENTRATION_NG_UL, CONCENTRATION_UNIT, VOLUME_UNIT, authorid)

    headers = {'Content-Type': 'application/json'}
    puppeteer_output_json = requests.post(PUPPETEER_URL, data=puppeteer_input_json, headers=headers)

    return tecan_json_to_gwl.get_tecan_instructions(puppeteer_output_json)


@app.route("/generate-autoprotocol", methods=['POST'])
def buildRequest():
    protocol_dict = tecan_to_autoprotocol.convert_tecan(request.get_json()) # returns instructions in dictionary form
    return json.dumps(protocol_dict)



if __name__ == "__main__":
    app.run(debug=True)
