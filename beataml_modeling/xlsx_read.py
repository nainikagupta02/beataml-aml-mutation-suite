"""Read XLSX cell values with Python's standard library, without editing workbooks."""
import zipfile, xml.etree.ElementTree as ET, posixpath, re
import pandas as pd
import numpy as np
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
NS={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
def workbook(name):
    out={}
    with zipfile.ZipFile(ROOT/name) as z:
        assert z.testzip() is None
        ss=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            ss=[''.join(s.itertext()) for s in ET.fromstring(z.read('xl/sharedStrings.xml'))]
        rel={r.attrib['Id']:r.attrib['Target'] for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        for sh in ET.fromstring(z.read('xl/workbook.xml')).find('m:sheets',NS):
            target=rel[sh.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']]
            path=posixpath.normpath('xl/'+target) if not target.startswith('/') else target[1:]
            rows=[]
            for row in ET.fromstring(z.read(path)).findall('m:sheetData/m:row',NS):
                vals={}
                for c in row:
                    col=re.match('[A-Z]+',c.attrib['r'])[0]; idx=0
                    for ch in col: idx=idx*26+ord(ch)-64
                    v=c.find('m:v',NS); typ=c.get('t'); val=None
                    if typ=='inlineStr': val=''.join(c.find('m:is',NS).itertext())
                    elif v is not None:
                        val=v.text
                        if typ=='s': val=ss[int(val)]
                        elif typ not in ('str','e'): val=float(val) if '.' in val or 'E' in val.upper() else int(val)
                    vals[idx-1]=val
                if any(v is not None for v in vals.values()): rows.append(vals)
            width=max(max(r) for r in rows)+1
            arr=[[r.get(i) for i in range(width)] for r in rows]
            out[sh.get('name')]=pd.DataFrame(arr[1:],columns=arr[0])
    return out
