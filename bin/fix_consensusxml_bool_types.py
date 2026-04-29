#!/usr/bin/env python
"""
Fix consensusXML files with invalid 'bool' UserParam types.
OpenMS only supports int type with 0/1 values, not 'bool'.
"""

import sys
import re

def fix_consensusxml_bool_types(input_file, output_file=None):
    """Convert type='bool' to type='int' with 0/1 values."""
    
    if output_file is None:
        output_file = input_file
    
    print(f"Reading: {input_file}")
    
    with open(input_file, 'r') as f:
        content = f.read()
    
    # Pattern to match: <UserParam type="bool" name="..." value="True|False"/>
    pattern = r'<UserParam type="bool" name="([^"]*)" value="([^"]*)"/>'
    
    def replacement(match):
        name = match.group(1)
        value = match.group(2)
        # Convert True/False to 1/0
        int_value = "1" if value.lower() == "true" else "0"
        return f'<UserParam type="int" name="{name}" value="{int_value}"/>'
    
    # Replace all bool types
    new_content = re.sub(pattern, replacement, content)
    
    # Count replacements
    num_replacements = len(re.findall(pattern, content))
    
    # Write output
    with open(output_file, 'w') as f:
        f.write(new_content)
    
    print(f"✓ Fixed {num_replacements} bool UserParam entries")
    print(f"✓ Saved to: {output_file}")
    
    return num_replacements

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python fix_consensusxml_bool_types.py <consensusxml_file> [output_file]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else input_file
    
    try:
        num_fixed = fix_consensusxml_bool_types(input_file, output_file)
        if num_fixed == 0:
            print("⚠ No bool types found - file may already be fixed")
    except Exception as e:
        print(f"✗ Error: {e}")
        sys.exit(1)
