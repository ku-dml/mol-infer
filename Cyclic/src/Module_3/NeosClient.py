# Copyright (c) 2017 NEOS-Server
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
Python source code - Python XML-RPC client for NEOS Server
"""

import argparse
import os
import sys
import time
try:
  import xmlrpc.client as xmlrpclib
except ImportError:
  import xmlrpclib
from collections import defaultdict

class NeosResult(object):
    def __init__(self, v):
        self._value = v

    def value(self):
        return self._value

    def __getitem__(self, item):
        return self


def make_xml(email, lp_string):
    return ('<document>\n'
            '<category>lp</category>\n'
            '<solver>CPLEX</solver>\n'
            '<inputMethod>LP</inputMethod>\n'
            '<email>' + email + '</email>\n'
            '<LP><![CDATA[\n' + lp_string + ']]></LP>\n'
            '<post><![CDATA[\ndisplay solution variables -\n]]></post>\n'
            '</document>')

def solveNeosMILP(email, lp_file):
    server = "https://neos-server.org:3333"

    neos = xmlrpclib.ServerProxy(server)

    alive = neos.ping()
    if alive != "NeosServer is alive\n":
        sys.stderr.write("Could not make connection to NEOS Server\n")
        sys.exit(1)

    with open(lp_file, 'r') as f:
        lp_string = f.read()
    xml = make_xml(email, lp_string)

    (jobNumber, password) = neos.submitJob(xml)
    sys.stdout.write("Job number = %d\nJob password = %s\n" % (jobNumber, password))
    if jobNumber == 0:
        sys.stderr.write("NEOS Server error: %s\n" % password)
        sys.exit(1)
    else:
        offset = 0
        status = ""
        while status != "Done":
            time.sleep(1)
            #(msg, offset) = neos.getIntermediateResults(jobNumber, password, offset)
            #sys.stdout.write(msg.data.decode())
            status = neos.getJobStatus(jobNumber, password)
        msg = neos.getFinalResults(jobNumber, password)
        result_str = msg.data.decode()
        #sys.stdout.write(result_str)

        # convert result string to dict
        variables_dict = defaultdict(lambda: NeosResult(0))
        skip = True
        for line in iter(result_str.splitlines()):
            if skip:    # skip solver logs
                if line.startswith('Variable Name'):
                    skip = False
                continue
            # parse line
            ls = line.split()
            if len(ls) > 2: # end of variables
                break
            _name_str, _value_str = ls
            if _name_str.startswith('ann_'):
                _name_str = _name_str[4:-4]
            parentheses_position = _name_str.find('(')  # first parentheses
            if parentheses_position == -1:  # if no parentheses found
                variables_dict[_name_str] = NeosResult(eval(_value_str))
            else:
                # separate name and value from line
                _name = _name_str[:parentheses_position]
                _key_str = _name_str[parentheses_position+1:]
                # do some cleaning
                _key_str = _key_str.replace(')(',',')
                _key_str = _key_str.replace('_','')
                # sometimes the last ')' is missing
                if _key_str.count('(') < _key_str.count(')'):
                    _key_str = _key_str[:-1]
                # variable names may contain different numbers of parentheses
                while _key_str.startswith('(') and _key_str.endswith(')'):
                    _key_str = _key_str[1:-1]
                _keys = _key_str.split(',')
                if len(_keys) == 1:
                    try:
                        _key = eval(_keys[0])
                    except NameError:
                        _key = eval('\''+_keys[0]+'\'')
                    except:
                        raise
                else:
                    _key_list = []
                    for _k in _keys:
                        # trim and count parentheses
                        lead_parentheses = 0
                        tail_parentheses = 0
                        while _k.startswith('('):
                            _k = _k[1:]
                            lead_parentheses += 1
                        while _k.endswith(')'):
                            _k = _k[:-1]
                            tail_parentheses += 1
                        # check if quotation marks are needed
                        try:
                            eval(_k)
                        except NameError:
                            _k = '\''+_k+'\''
                        except:
                            raise
                        # add parentheses back
                        _k = '('*lead_parentheses + _k + ')'*tail_parentheses
                        _key_list.append(_k)
                    _key = eval('(' + ','.join(_key_list) + ')')
                #print(_name, type(_key), _key, eval(_value_str))
                if _name not in variables_dict:
                    variables_dict[_name] = defaultdict(lambda: NeosResult(0))
                variables_dict[_name][_key] = NeosResult(eval(_value_str))
    return variables_dict

    # vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

