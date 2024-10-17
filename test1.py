class Elliptic_curves:
    def __init__(self, a, b, p):
        self.a = a
        self.b = b
        self.p = p

    def is_on_curve(self,x,y):
        """Check if the point is on the curve."""
        left = (y ** 2) % self.p
        right = (x ** 3 + self.a * x + self.b) % self.p
        
        assert left == right, f"point is not on the curve" 

def point_addition(P, Q):
    """Add two points P and Q on the elliptic curve."""
    
    x1,y1 = P
    x2,y2 = Q
    
    if P == Q:  # Point doubling
        return point_doubling(P)
    
    assert x1 != x2 and y1 != -y2 % curev.p, f"This addition will return point at infinity"
    
    # Compute the slope
    m = ((y2 - y1) * manual_pow((x2 - x1), -1, curev.p)) % curev.p

    # Compute the resulting point
    x_r = (m ** 2 - x2 - x1) % curev.p
    y_r = (m * (x1 - x_r) - y1) % curev.p

    return (x_r, y_r)

def point_doubling(P):
    """Double a point P on the elliptic curve."""
    x1,y1 = P
    m = ((3 * x1**2 + curev.a) * manual_pow(2 * y1, -1, curev.p)) % curev.p
    x_r = (m**2 - 2 * x1) % curev.p
    y_r = (m * (x1 - x_r) - y1) % curev.p
    return (x_r, y_r)

def scalar_multiplication(P, k):
    """Multiply a point P by an integer k."""
    result = None
    addend = P

    while k:
        if k & 1:
            if result is None:
                result = addend
            else:
                result = point_addition(result, addend)

        addend = point_doubling(addend)
        k >>= 1 #shift right

    return result

def manual_pow(base, exp, mod):
    # If exp is negative, calculate the modular inverse of the base
    if exp < 0:
        # Calculate modular inverse of base using Fermat's Little Theorem
        base = manual_pow(base, mod - 2, mod)
        exp = -exp  # Convert exponent to positive

    result = 1  # Initialize result to 1
    base = base % mod  # Take mod of base to handle cases where base >= mod
    
    while exp > 0:
        # If exp is odd, multiply base with result and take mod
        if exp % 2 == 1:
            result = (result * base) % mod
        
        # Now exp is even, divide it by 2
        exp = exp // 2
        #print(result)
        
        # Square the base and take mod
        base = (base * base) % mod
    return result





selector =int(input("Enter the curve selector"))
if( selector == 1 ) :
    a    = 2
    b    = 2
    p    = 17
    G1_x = 5
    G1_y = 1
    G1 = (G1_x, G1_y) 
else : 
    a    = 0
    b    = 1 
    p    = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eab1fffeb153ffffb9feffffffffaaab
    G1_x = 0x17F1D3A73197D7942695638C4FA9AC0FC3688C4F9774B905A14E3A3F171BAC586C55E83FF97A1AEFFB3AF00ADB22C6BB
    G1_y = 0x08B3F481E3AAA0F1A09E30ED741D8AE4FCF5E095D5D00AF600DB18CB2C04B3EDD03CC744A2888AE40CAA232946C5E7E1
    G1 = (G1_x, G1_y)   

curev = Elliptic_curves(a,b,p)

print("G1 on curve:", curev.is_on_curve(G1_x,G1_y))

# Perform scalar multiplication
k = 6
result = scalar_multiplication(G1, k)
print(result)
#print(f"Result of {k} * G1: ({result.x}, {result.y})")