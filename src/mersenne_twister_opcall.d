/++
The Mersenne Twister generator.

Copyright: Copyright Andrei Alexandrescu 2008 - 2009, Ilya Yaroshenko 2016-.
License:   $(HTTP www.boost.org/LICENSE_1_0.txt, Boost License 1.0).
Authors: $(HTTP erdani.org, Andrei Alexandrescu) Ilya Yaroshenko (rework)
+/
module mersenne_twister_opcall;

// Segments of the code in this file Copyright (c) 1997 by Rick Booth
// From "Inner Loops" by Rick Booth, Addison-Wesley

// Work derived from:

/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

import std.traits;

/++
The $(LUCKY Mersenne Twister) generator.
+/
struct MersenneTwisterEngine(UIntType, size_t w, size_t n, size_t m, size_t r,
                             UIntType a, size_t u, UIntType d, size_t s,
                             UIntType b, size_t t,
                             UIntType c, size_t l, UIntType f)
    if (isUnsigned!UIntType)
{
    ///
    enum isRandomEngine = true;

    static assert(0 < w && w <= UIntType.sizeof * 8);
    static assert(1 <= m && m <= n);
    static assert(0 <= r && 0 <= u && 0 <= s && 0 <= t && 0 <= l);
    static assert(r <= w && u <= w && s <= w && t <= w && l <= w);
    static assert(0 <= a && 0 <= b && 0 <= c);

    @disable this();
    @disable this(this);

    /// Largest generated value.
    enum UIntType max = UIntType.max >> (UIntType.sizeof * 8u - w);
    static assert(a <= max && b <= max && c <= max && f <= max);

    private enum UIntType lowerMask = (cast(UIntType) 1u << r) - 1;
    private enum UIntType upperMask = ~lowerMask & max;

    /**
    Parameters for the generator.
    */
    enum size_t   wordSize   = w;
    enum size_t   stateSize  = n; /// ditto
    enum size_t   shiftSize  = m; /// ditto
    enum size_t   maskBits   = r; /// ditto
    enum UIntType xorMask    = a; /// ditto
    enum size_t   temperingU = u; /// ditto
    enum UIntType temperingD = d; /// ditto
    enum size_t   temperingS = s; /// ditto
    enum UIntType temperingB = b; /// ditto
    enum size_t   temperingT = t; /// ditto
    enum UIntType temperingC = c; /// ditto
    enum size_t   temperingL = l; /// ditto
    enum UIntType initializationMultiplier = f; /// ditto


    /// The default seed value.
    enum UIntType defaultSeed = 5489;

    /++
    Current reversed payload index with initial value equals to `n-1`
    +/
    size_t index = void;

    private UIntType _z = void;

    /++
    Reversed(!) payload.
    +/
    UIntType[n] data = void;

    /++
    Constructs a MersenneTwisterEngine object.
    +/
    this(UIntType value) @safe pure nothrow @nogc
    {
        static if (max == UIntType.max)
            data[$-1] = value;
        else
            data[$-1] = value & max;
        foreach_reverse (size_t i, ref e; data[0 .. $-1])
        {
            e = f * (data[i + 1] ^ (data[i + 1] >> (w - 2))) + cast(UIntType)(n - (i + 1));
            static if (max != UIntType.max)
                e &= max;
        }
        index = n-1;
        opCall();
    }

    /++
    Advances the generator.
    +/
    UIntType opCall() @safe pure nothrow @nogc
    {
        // This function blends two nominally independent
        // processes: (i) calculation of the next random
        // variate from the cached previous `data` entry
        // `_z`, and (ii) updating `data[index]` and `_z`
        // and advancing the `index` value to the next in
        // sequence.
        //
        // By interweaving the steps involved in these
        // procedures, rather than performing each of
        // them separately in sequence, the variables
        // are kept 'hot' in CPU registers, allowing
        // for significantly faster performance.
        sizediff_t index = this.index;
        sizediff_t next = index - 1;
        if(next < 0)
            next = n - 1;
        auto z = _z;
        sizediff_t conj = index - m;
        if(conj < 0)
            conj = index - m + n;
        static if (d == UIntType.max)
            z ^= (z >> u);
        else
            z ^= (z >> u) & d;
        auto q = data[index] & upperMask;
        auto p = data[next] & lowerMask;
        z ^= (z << s) & b;
        auto y = q | p;
        auto x = y >> 1;
        z ^= (z << t) & c;
        if (y & 1)
            x ^= a;
        auto e = data[conj] ^ x;
        z ^= (z >> l);
        _z = data[index] = e;
        this.index = next;
        return z;
    }
}

/++
A $(D MersenneTwisterEngine) instantiated with the parameters of the
original engine $(HTTP en.wikipedia.org/wiki/Mersenne_Twister,
MT19937), generating uniformly-distributed 32-bit numbers with a
period of 2 to the power of 19937.

This is recommended for random number generation on 32-bit systems
unless memory is severely restricted, in which case a
$(REF_ALTTEXT Xorshift, Xorshift, mir, random, engine, xorshift)
would be the generator of choice.
+/
alias Mt19937 = MersenneTwisterEngine!(uint, 32, 624, 397, 31,
                                       0x9908b0df, 11, 0xffffffff, 7,
                                       0x9d2c5680, 15,
                                       0xefc60000, 18, 1812433253);

///
@safe unittest
{
    import mir.random.engine;

    // bit-masking by generator maximum is necessary
    // to handle 64-bit `unpredictableSeed`
    auto gen = Mt19937(unpredictableSeed & Mt19937.max);
    auto n = gen();

    import std.traits;
    static assert(is(ReturnType!gen == uint));
}

/++
A $(D MersenneTwisterEngine) instantiated with the parameters of the
original engine $(HTTP en.wikipedia.org/wiki/Mersenne_Twister,
MT19937), generating uniformly-distributed 64-bit numbers with a
period of 2 to the power of 19937.

This is recommended for random number generation on 64-bit systems
unless memory is severely restricted, in which case a
$(REF_ALTTEXT Xorshift, Xorshift, mir, random, engine, xorshift)
would be the generator of choice.
+/
alias Mt19937_64 = MersenneTwisterEngine!(ulong, 64, 312, 156, 31,
                                          0xb5026f5aa96619e9, 29, 0x5555555555555555, 17,
                                          0x71d67fffeda60000, 37,
                                          0xfff7eee000000000, 43, 6364136223846793005);
