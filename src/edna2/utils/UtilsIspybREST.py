#
# Copyright (c) European Synchrotron Radiation Facility (ESRF)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

__authors__ = ["A. Finke"]
__license__ = "MIT"
__date__ = "11/10/2023"

import os
import json
import time
import requests
from datetime import datetime
from pathlib import Path

import requests

from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging

logger = UtilsLogging.getLogger()


def authenticate():
    """returns ISPyB Token"""
    config = UtilsConfig.getTaskConfig("ISPyB")
    proxies = {
        "http_proxy": config.get("http_proxy"),
        "https_proxy": config.get("https_proxy"),
    }
    site = config.get("site")
    url = config.get("rest_url", "http://localhost") + "/authenticate?site=" + site
    try:
        r = requests.post(
            url,
            headers={"content-type": "application/x-www-form-urlencoded"},
            proxies=proxies,
            data={"login": config.get("username"), "password": config.get("password")},
            verify=True
        )
        token = json.loads(r.text)["token"]
    except:
        token = None
    return token

def getPhasingViewByDataCollectionId(proposal, dataCollectionId, token=None):
    if token is None:
        token = authenticate()

    config = UtilsConfig.getTaskConfig("ISPyB")
    root_url = config.get("rest_url")

    phasing = f"{token}/proposal/MX{proposal}/mx/phasing/datacollectionid/{dataCollectionId}/list"
    r = requests.get(root_url + phasing)
    try:
        output = json.loads(r.text)
        output = output[0][0]
    except:
        output = None
    return output
