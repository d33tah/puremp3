import sys
from bitarray import bitarray
from bitarray.util import int2ba

def main():
    # Initialize a bitarray to store the reconstructed bits
    bit_writer = bitarray()

    # Read input lines
    input_lines = sys.stdin.read().strip().split("\n")

    for line in input_lines:
        if line.startswith("read_bit()"):
            # Extract the bit (assume 'true' for 1 and 'false' for 0)
            value = line.split(": ")[1].strip().lower() == "true"
            bit_writer.append(value)
        
        elif line.startswith("read_bits("):
            # Extract number of bits to read and the value
            n, value = line.split("): ")
            num_bits = int(n.split("(")[1])
            value = int(value.strip())
            # Append the corresponding bitarray
            bits = int2ba(value, length=num_bits, endian='big')
            bit_writer.extend(bits)
        
        elif line.startswith("skip_bits("):
            # Extract number of bits to skip
            n = line.split("(")[1].split(")")[0]
            num_bits = int(n.strip())
            # Extend with zeros for skipped bits
            bit_writer.extend([False] * num_bits)

    # Write the bitarray to a binary file
    with open("a.bin", "wb") as f:
        bit_writer.tofile(f)

if __name__ == "__main__":
    main()
