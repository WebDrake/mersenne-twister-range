/++
The Mersenne Twister generator.

Copyright: Copyright Andrei Alexandrescu 2008 - 2009, Ilya Yaroshenko 2016-.
License:   $(HTTP www.boost.org/LICENSE_1_0.txt, Boost License 1.0).
Authors: $(HTTP erdani.org, Andrei Alexandrescu) Ilya Yaroshenko (rework)
+/
module mersenne_twister_range;

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
                             UIntType a,
                             uint u, UIntType d,
                             uint s, UIntType b,
                             uint t, UIntType c,
                             uint l)
    if (isUnsigned!UIntType)
{
    import std.range.primitives;

    static assert(0 < w && w <= UIntType.sizeof * 8);
    static assert(1 <= m && m <= n);
    static assert(0 <= r && 0 <= u && 0 <= s && 0 <= t && 0 <= l);
    static assert(r <= w && u <= w && s <= w && t <= w && l <= w);
    static assert(0 <= a && 0 <= b && 0 <= c);
    static assert(n < UIntType.max);

    enum bool isUniformRandom = true;

    private enum UIntType upperMask = ~((cast(UIntType) 1u << (UIntType.sizeof * 8 - (w - r))) - 1);
    private enum UIntType lowerMask = (cast(UIntType) 1u << r) - 1;

    /**
    Parameters for the generator.
    */
    enum size_t wordSize   = w;
    enum size_t stateSize  = n; /// ditto
    enum size_t shiftSize  = m; /// ditto
    enum size_t maskBits   = r; /// ditto
    enum UIntType xorMask    = a; /// ditto
    enum uint temperingU     = u; /// ditto
    enum UIntType temperingD = d; /// ditto
    enum uint temperingS     = s; /// ditto
    enum UIntType temperingB = b; /// ditto
    enum uint temperingT     = t; /// ditto
    enum UIntType temperingC = c; /// ditto
    enum uint temperingL     = l; /// ditto

    /// Smallest generated value (0)
    enum UIntType min = 0;

    /// Largest generated value.
    enum UIntType max = UIntType.max >> (UIntType.sizeof * 8u - w);
    static assert(a <= max && b <= max && c <= max);

    /// The default seed value.
    enum UIntType defaultSeed = 5489;

    /// Collection of all state variables
    /// used by the generator
    private struct State
    {
        private UIntType y = void;
        private UIntType z = 0;
        private UIntType index = void;
        private UIntType[n] data = void;

    }

    /// State variables used by the generator;
    /// initialized to values equivalent to
    /// explicitly seeding the generator with
    /// `defaultSeed`
    private State state = defaultState();
    // NOTE: the above is a workaround to ensure
    // backwards compatibility with the original
    // implementation, which permitted implicit
    // construction.  With `@disable this();`
    // it would not be necessary.

    /**
       Constructs a MersenneTwisterEngine object.
    */
    this(UIntType value) @safe pure nothrow @nogc
    {
        this.seed(value);
    }

    /**
       Constructs a MersenneTwisterEngine object using a range
       to seed the generator, which must have at least as many
       elements as `data`.
    */
    this(T)(T range)
        if (isInputRange!T && is(Unqual!(ElementType!T) == UIntType))
    {
        this.seed(range);
    }

    /**
       Constructs a MersenneTwisterEngine object by
       duplicating the state of an existing instance
    */
    this(ref typeof(this) that)
    {
        this.state = that.state;
    }

    /**
       Generates the default initial state for a Mersenne
       Twister; equivalent to the internal state obtained
       by calling `seed(defaultSeed)`
    */
    private static State defaultState() @safe pure nothrow @nogc
    {
        State mtState;
        seedImpl(defaultSeed, mtState);
        return mtState;
    }

    /**
       (Re)seeds the generator
    */
    void seed(UIntType value) @safe pure nothrow @nogc
    {
        this.seedImpl(value, this.state);
    }

    /**
       Implementation of the seeding mechanism, which
       can be used with an arbitrary `State` instance
    */
    private static void seedImpl(UIntType value, ref State mtState)
    {
        static if (w == UIntType.sizeof * 8)
        {
            mtState.data[$-1] = value;
        }
        else
        {
            static assert(max + 1 > 0);
            mtState.data[$-1] = value % (max + 1);
        }
        static if (is(UIntType == uint))
            enum UIntType f = 1812433253;
        else static if (is(UIntType == ulong))
            enum UIntType f = 6364136223846793005;
        else
            static assert(0, UIntType.stringof ~ " is not supported by MersenneTwisterEngine.");
        foreach_reverse (size_t i, ref e; mtState.data[0 .. $-1])
            e = f * (mtState.data[i + 1] ^ (mtState.data[i + 1] >> (w - 2))) + cast(UIntType)(n - (i + 1));
        mtState.index = n-1;

        // double popFront() to guarantee both
        // `_z` and `_y` are derived from the
        // newly set values in `data`
        MersenneTwisterEngine.popFrontImpl(mtState);
        MersenneTwisterEngine.popFrontImpl(mtState);
    }

    /**
       (Re)seeds the generator using values from a range,
       which must have at least as many elements as `data`
    */
    void seed(T)(T range)
        if (isInputRange!T && is(Unqual!(ElementType!T) == UIntType))
    {
        this.seedImpl(range, this.state);
    }

    /**
       Implementation of the range-based seeding mechanism,
       which can be used with an arbitrary `State` instance
    */
    private static void seedImpl(T)(T range, ref State mtState)
        if (isInputRange!T && is(Unqual!(ElementType!T) == UIntType))
    {
        size_t j;
        for (j = 0; j < n && !range.empty; ++j, range.popFront())
        {
            sizediff_t idx = n - j - 1;
            mtState.data[idx] = range.front;
        }

        mtState.index = n - 1;

        if (range.empty && j < n)
        {
            import core.internal.string : UnsignedStringBuf, unsignedToTempString;

            UnsignedStringBuf buf = void;
            string s = "MersenneTwisterEngine.seed: Input range didn't provide enough elements: Need ";
            s ~= unsignedToTempString(n, buf, 10) ~ " elements.";
            throw new Exception(s);

        }

        // double popFront() to guarantee both
        // `_z` and `_y` are derived from the
        // newly set values in `data`
        MersenneTwisterEngine.popFrontImpl(mtState);
        MersenneTwisterEngine.popFrontImpl(mtState);
    }

    /// Range primitive: a Mersenne Twister is never empty
    enum bool empty = false;

    /// Range primitive: get the current random variate
    UIntType front() @property @safe pure nothrow @nogc
    {
        return this.state.y;
    }

    /// Range primitive: advance the generator state
    /// to get the next random variate
    void popFront() @safe pure nothrow @nogc
    {
        this.popFrontImpl(this.state);
    }

    /// Internal implementation of `popFront()`, which
    /// can be used with an arbitrary `State` instance
    private static void popFrontImpl(ref State mtState)
    {
        // This function blends two nominally independent
        // processes: (i) calculation of the next random
        // variate `mtState.y` from the cached previous
        // `data` entry `z`, and (ii) updating the value
        // of `data[index]` and `mtState.z` and advancing
        // the `index` value to the next in sequence.
        //
        // By interweaving the steps involved in these
        // procedures, rather than performing each of
        // them separately in sequence, the variables
        // are kept 'hot' in CPU registers, allowing
        // for significantly faster performance.
        sizediff_t index = cast(size_t)mtState.index;
        sizediff_t next = index - 1;
        if (next < 0)
            next = n - 1;
        auto z = mtState.z;
        sizediff_t conj = index - m;
        if (conj < 0)
            conj = index - m + n;
        static if (d == UIntType.max)
            z ^= (z >> u);
        else
            z ^= (z >> u) & d;
        auto q = mtState.data[index] & upperMask;
        auto p = mtState.data[next] & lowerMask;
        z ^= (z << s) & b;
        auto y = q | p;
        auto x = y >> 1;
        z ^= (z << t) & c;
        if (y & 1)
            x ^= a;
        auto e = mtState.data[conj] ^ x;
        z ^= (z >> l);
        mtState.z = mtState.data[index] = e;
        mtState.index = cast(UIntType)next;
        mtState.y = z;
    }

    /// Range primitive: return an identical copy
    /// of the generator
    typeof(this) save() @property @safe pure nothrow @nogc
    {
        return typeof(this)(this);
    }
}

/++
A $(D MersenneTwisterEngine) instantiated with the parameters of the
original engine $(HTTP en.wikipedia.org/wiki/Mersenne_Twister,
MT19937), generating uniformly-distributed 32-bit numbers with a
period of 2 to the power of 19937.
+/
alias Mt19937_32 = MersenneTwisterEngine!(uint, 32, 624, 397, 31,
                                       0x9908b0df,
                                       11, 0xffffffff,
                                        7, 0x9d2c5680,
                                       15, 0xefc60000,
                                       18);
/++
A $(D MersenneTwisterEngine) instantiated with the parameters of the
original engine $(HTTP en.wikipedia.org/wiki/Mersenne_Twister,
MT19937), generating uniformly-distributed 64-bit numbers with a
period of 2 to the power of 19937.
+/
alias Mt19937_64 = MersenneTwisterEngine!(ulong, 64, 312, 156, 31,
                                       0xb5026f5aa96619e9,
                                       29, 0x5555555555555555,
                                       17, 0x71d67fffeda60000,
                                       37, 0xfff7eee000000000,
                                       43);

unittest
{
    import std.random : isUniformRNG;

    static assert(isUniformRNG!Mt19937_32);
    static assert(isUniformRNG!Mt19937_64);
}

unittest
{
    import std.random : Mt19937;
    import std.range : iota, take;

    auto gen = Mt19937_32(iota(1000u));

    Mt19937 gen_std;
    gen_std.seed(iota(1000u));

    foreach (_; 0 .. 1000)
    {
        assert(gen.front == gen_std.front);
        gen.popFront();
        gen_std.popFront();
    }
}

unittest
{
    auto genA = Mt19937_32(101);
    auto genB = genA.save;

    foreach (_; 0 .. 1000)
    {
        assert(genA.front == genB.front);
        genA.popFront();
        genB.popFront();
    }
}

unittest
{
    import std.random : Mt19937;
    import std.range : iota, popFrontN, take;
    import std.stdio : writeln;

    {
        Mt19937 gen_std;
        gen_std.take(10).writeln;
    }

    {
        Mt19937_32 gen;
        (&gen).take(10).writeln;
    }

    {
        Mt19937_32 gen = Mt19937_32(Mt19937_32.defaultSeed);
        (&gen).take(10).writeln;
    }

    {
        Mt19937 gen_std;
        popFrontN(gen_std, 9999);
        writeln(gen_std.front);
    }

    {
        Mt19937_32 gen;
        popFrontN(gen, 9999);
        writeln(gen.front);
    }

    {
        Mt19937 gen_std;
        gen_std.seed(iota(1u, 625u));
        gen_std.take(10).writeln;
    }

    {
        Mt19937_32 gen;
        gen.seed(iota(1u, 625u));
        (&gen).take(10).writeln;
    }

    {
        Mt19937_32 genA;
        auto genB = Mt19937_32(Mt19937_32.defaultSeed);
        assert(genA.state == genB.state);
    }
}
