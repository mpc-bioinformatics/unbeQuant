#!/usr/bin/env python3
"""
Fix consensusXML files with invalid 'bool' UserParam types using proper XML parsing.
OpenMS only supports int type with 0/1 values, not 'bool'.
"""

import sys
import xml.etree.ElementTree as ET

def fix_consensusxml_bool_types_proper(input_file, output_file=None):
    """Convert type='bool' to type='int' with 0/1 values using proper XML parsing."""
    
    if output_file is None:
        output_file = input_file
    
    print(f"Reading: {input_file}")
    
    # Register namespaces to preserve them
    ET.register_namespace('xsi', 'http://www.w3.org/2001/XMLSchema-instance')
    
    # Parse the XML
    tree = ET.parse(input_file)
    root = tree.getroot()
    
    # Find all elements with type="bool"
    num_fixed = 0
    
    # Recursively search for all UserParam elements
    for elem in root.iter():
        if elem.tag == 'UserParam':
            if elem.get('type') == 'bool':
                # Convert bool to int
                elem.set('type', 'int')
                
                # Convert True/False to 1/0
                value = elem.get('value', '')
                if value.lower() == 'true':
                    elem.set('value', '1')
                elif value.lower() == 'false':
                    elem.set('value', '0')
                
                num_fixed += 1
    
    # Write output with proper XML declaration
    tree.write(output_file, encoding='ISO-8859-1', xml_declaration=True)
    
    print(f"✓ Fixed {num_fixed} bool UserParam entries")
    print(f"✓ Saved to: {output_file}")
    
    # Add back the XML stylesheet if it was there
    try:
        with open(output_file, 'r') as f:
            content = f.read()
        
        # Check if we need to add the stylesheet
        if 'xml-stylesheet' not in content:
            # Insert stylesheet after XML declaration
            lines = content.split('\n', 1)
            if len(lines) == 2:
                content = lines[0] + '\n<?xml-stylesheet type="text/xsl" href="https://www.openms.de/xml-stylesheet/ConsensusXML.xsl" ?>\n' + lines[1]
                with open(output_file, 'w') as f:
                    f.write(content)
                print("✓ Added XML stylesheet declaration")
    except Exception as e:
        print(f"⚠ Warning: Could not add stylesheet: {e}")
    
    return num_fixed

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 fix_consensusxml_bool_types_proper.py <consensusxml_file> [output_file]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else input_file
    
    try:
        num_fixed = fix_consensusxml_bool_types_proper(input_file, output_file)
        if num_fixed == 0:
            print("⚠ No bool types found - file may already be fixed")
    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
