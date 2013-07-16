import xml.etree.ElementTree as et

def read_xml(filename):
    
    def recursive_dict(element):
        return element.tag, \
            dict(map(recursive_dict, element)) or element.text
    
    
    doc = et.parse(filename)
    
    return recursive_dict(doc.getroot())
