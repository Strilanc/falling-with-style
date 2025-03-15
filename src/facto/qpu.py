class Qubit:
    """A qubit identifier."""
    def __init__(self, *, index: int):
        self.index = index
        self.debug_val = 0


class QPU:
    """A simple quantum circuit builder.

    Has basic support for doing reversible circuit sims to verify correctness of methods.
    """
    def __init__(self):
        self.x_count = 0
        self.cx_count = 0
        self.ccx_count = 0
        self.m_count = 0
        self.recorded = None
        self.next_qubit_index = 0
        self.allocated = []

    def alloc_qubit(self) -> Qubit:
        """Allocates a clean qubit."""
        if self.allocated:
            return self.allocated.pop()
        q = Qubit(index=self.next_qubit_index)
        self.next_qubit_index += 1
        return q

    def free_qubit(self, q: Qubit, *, may_not_be_zero: bool = False):
        """Deallocates a qubit."""
        if not may_not_be_zero:
            assert q.debug_val == 0
        q.debug_val = 0
        self.allocated.append(q)

    def debug_set_val(self, vals: list[Qubit], mi: int, val: int):
        """Debugging method for forcing tracked values.

        Args:
            vals: The qubit register to write to.
            mi: The instance index to write to.
            val: The value to write.
        """
        for k, v in enumerate(vals):
            v.debug_val &= ~(1 << mi)
            if val & (1 << k):
                v.debug_val ^= (1 << mi)

    def debug_read_val(self, vals: list[Qubit], mi: int) -> int:
        """Debugging method for reading tracked values.

        Args:
            vals: The qubit register to read from.
            mi: The instance index to read from.

        Returns:
            The read value.
        """
        total = 0
        for k, v in enumerate(vals):
            if v.debug_val & (1 << mi):
                total |= 1 << k
        return total

    def debug_read_vals(self, vals: list[Qubit], m: int) -> list[int]:
        """Debugging method for reading tracked values.

        Args:
            vals: The qubit register to read from.
            m: The number of instances to read up to.

        Returns:
            The read values.
        """
        return [self.debug_read_val(vals, k) for k in range(m)]

    def x(self, *, t: Qubit):
        """Bit flips a qubit."""
        t.debug_val ^= -1
        self.x_count += 1
        if self.recorded is not None:
            self.recorded.append(('X', t.index))

    def h(self, *, t: Qubit):
        """Applies Hadamard to a qubit."""
        if self.recorded is not None:
            self.recorded.append(('H', t.index))

    def m(self, *, t: Qubit):
        """Measures a qubit."""
        if self.recorded is not None:
            self.recorded.append(('M', t.index))
        self.m_count += 1

    def cphase(self, *, c1: Qubit, c2: Qubit, exponent: float):
        """Phases the 11 state of c1 and c2."""
        if self.recorded is not None:
            self.recorded.append(('CPHASE', exponent, c1.index, c2.index))

    def cx(self, *, c: Qubit | bool, t: Qubit):
        """Performs `t ^= c`."""
        if isinstance(c, bool):
            if c:
                self.x(t=t)
        else:
            t.debug_val ^= c.debug_val
            self.cx_count += 1
            if self.recorded is not None:
                self.recorded.append(('CX', c.index, t.index))

    def ccx(self, c1: Qubit, c2: Qubit, *, t: Qubit):
        """Performs `t ^= c1 & c2`."""
        t.debug_val ^= c1.debug_val & c2.debug_val
        self.ccx_count += 1
        if self.recorded is not None:
            self.recorded.append(('CCX', c1.index, c2.index, t.index))

    def iadd(self, *, target: list[Qubit], offset: list[Qubit]):
        """Performs `target += offset`.

        This is the Cuccaro adder from https://arxiv.org/pdf/quant-ph/0410184.
        """
        assert len(target) <= len(offset) + 1
        if len(target) < len(offset):
            offset = offset[:len(target)]
        n = len(offset)

        carry_in = self.alloc_qubit()
        offset = [carry_in, *offset]

        # MAJ sweep.
        for k in range(n):
            a = offset[k]
            b = target[k]
            c = offset[k + 1]
            self.cx(c=c, t=a)
            self.cx(c=c, t=b)
            self.ccx(c1=a, c2=b, t=c)

        # Carry out.
        if len(target) == len(offset):
            self.cx(t=target[-1], c=offset[-1])

        # UMA sweep.
        for k in range(n)[::-1]:
            a = offset[k]
            b = target[k]
            c = offset[k + 1]
            self.ccx(c1=a, c2=b, t=c)
            self.cx(c=c, t=a)
            self.cx(c=a, t=b)

        self.free_qubit(carry_in)

    def flip_if_lt(self, *, lhs: list[Qubit], rhs: list[Qubit], out: Qubit):
        """Performs `out ^= lhs < rhs`.

        Uses a tweaked Cuccaro adder.
        """
        assert len(lhs) == len(rhs), f"flip_if_lt {len(lhs)=} != {len(rhs)=}"
        n = len(rhs)

        carry_in = self.alloc_qubit()
        rhs = [carry_in, *rhs]
        for e in lhs:
            self.x(t=e)

        # MAJ sweep.
        for k in range(n):
            a = rhs[k]
            b = lhs[k]
            c = rhs[k + 1]
            self.cx(c=c, t=a)
            self.cx(c=c, t=b)
            self.ccx(c1=a, c2=b, t=c)

        # Carry out.
        self.cx(t=out, c=rhs[-1])

        # Un-MAJ sweep.
        for k in range(n)[::-1]:
            a = rhs[k]
            b = lhs[k]
            c = rhs[k + 1]
            self.ccx(c1=a, c2=b, t=c)
            self.cx(c=c, t=a)
            self.cx(c=c, t=b)
        for e in lhs:
            self.x(t=e)
        self.free_qubit(carry_in)

    def isub(self, *, target: list[Qubit], offset: list[Qubit]):
        """Performs `target -= offset`"""
        for t in target:
            self.x(t=t)
        self.iadd(target=target, offset=offset)
        for t in target:
            self.x(t=t)

    def iadd_const_if(self, *, target: list[Qubit], control: Qubit | bool, offset_const: int):
        """Performs `if control: target += offset_const`"""
        offset = [self.alloc_qubit() for _ in range(len(target))]
        for k in range(len(target)):
            if offset_const & (1 << k):
                self.cx(c=control, t=offset[k])
        self.iadd(target=target, offset=offset)
        for k in range(len(target)):
            if offset_const & (1 << k):
                self.cx(c=control, t=offset[k])
        for q in offset:
            self.free_qubit(q)

    def flip_if_lt_const(self, *, lhs: list[Qubit], rhs_const: int, out: Qubit, control: bool | Qubit = True):
        """Performs `if control: out ^= lhs < rhs_const`"""
        rhs = [self.alloc_qubit() for _ in range(len(lhs))]
        for k in range(len(lhs)):
            if rhs_const & (1 << k):
                self.cx(c=control, t=rhs[k])
        self.flip_if_lt(lhs=lhs, rhs=rhs, out=out)
        for k in range(len(lhs)):
            if rhs_const & (1 << k):
                self.cx(c=control, t=rhs[k])
        for q in rhs:
            self.free_qubit(q)

    def flip_if_ge_const(self, *, lhs: list[Qubit], rhs_const: int, out: Qubit, control: bool | Qubit = True):
        """Performs `if control: out ^= lhs >= rhs_const`"""
        self.flip_if_lt_const(lhs=lhs, rhs_const=rhs_const, out=out, control=control)
        self.cx(c=control, t=out)

    def isub_const(self, *, target: list[Qubit], control: Qubit, offset_const: int):
        """Performs `if control: target -= offset_cond`"""
        self.iadd_const_if(target=target, control=control, offset_const=-offset_const)

    def iadd_const_mod_if(self, *, target: list[Qubit], control: Qubit, offset_const: int, modulus: int):
        """Performs `if control: target += offset_cond (mod modulus)`"""
        offset_const %= modulus
        neg_offset = (-offset_const) % modulus
        cmp_helper = self.alloc_qubit()
        self.isub_const(target=target + [cmp_helper], control=control, offset_const=neg_offset)
        self.iadd_const_if(target=target, control=cmp_helper, offset_const=modulus)
        self.flip_if_ge_const(lhs=target, rhs_const=modulus - neg_offset, out=cmp_helper, control=control)
        self.free_qubit(cmp_helper)

    def imul_const_mod_if(
            self,
            *,
            target: list[Qubit],
            const_factor: int,
            modulus: int,
            control: Qubit,
    ):
        """Performs `if control: target *= const_factor (mod modulus)`"""
        helper = [self.alloc_qubit() for _ in target]

        # Compute the product out of place.
        combined_control = self.alloc_qubit()
        for k in range(len(target)):
            self.ccx(c1=control, c2=target[k], t=combined_control)
            self.iadd_const_mod_if(
                target=helper,
                control=combined_control,
                offset_const=const_factor << k,
                modulus=modulus,
            )
            self.ccx(c1=control, c2=target[k], t=combined_control)
        self.free_qubit(combined_control)

        # Uncompute the input.
        inv_const_factor = -pow(const_factor, -1, modulus)
        for k in range(len(target)):
            self.iadd_const_mod_if(
                target=target,
                control=helper[k],
                offset_const=inv_const_factor << k,
                modulus=modulus,
            )

        # Swap output back into the target register.
        for k in range(len(target)):
            self.cx(c=helper[k], t=target[k])
            self.ccx(t=helper[k], c1=target[k], c2=control)

        for q in helper:
            self.free_qubit(q)

    def inv_qft_measure(self, target: list[Qubit]):
        """Performs a frequency basis measurement."""
        for k in range(len(target))[::-1]:
            q = target[k]
            for k2 in range(k + 1, len(target)):
                self.cphase(c1=q, c2=target[k2], exponent=-2**-(k2 - k))
            self.h(t=q)
            self.m(t=q)

    def perform_shor_quantum_part(self, *, g: int, modulus: int):
        n = modulus.bit_length()

        # Initialize the phase estimation register.
        phase_reg = []
        for k in range(2 * n):
            q = self.alloc_qubit()
            self.h(t=q)
            phase_reg.append(q)

        # Compute the modular exponentation.
        work_reg = []
        for k in range(n):
            work_reg.append(self.alloc_qubit())
        self.x(t=work_reg[0])
        for k in range(len(phase_reg)):
            q = phase_reg[k]
            self.imul_const_mod_if(
                target=work_reg,
                const_factor=pow(g, 2**k, modulus),
                modulus=modulus,
                control=q,
            )

        # Measure the frequency.
        self.inv_qft_measure(phase_reg)

        for q in phase_reg:
            self.free_qubit(q, may_not_be_zero=True)
        for q in work_reg:
            self.free_qubit(q, may_not_be_zero=True)
        assert len(self.allocated) == self.next_qubit_index

    def to_qasm(self, *, cphase_angle_cutoff: float = 0) -> str:
        assert self.recorded is not None
        lines = []
        lines.append('OPENQASM 2.0;')
        lines.append('include "qelib1.inc";')
        lines.append(f'qreg q[{self.next_qubit_index}];')
        lines.append(f'creg c[{self.m_count}];')
        m = 0

        for e in self.recorded:
            if e[0] == 'X':
                _, a = e
                lines.append(f'x q[{a}];')
            elif e[0] == 'CX':
                _, a, b = e
                lines.append(f'cx q[{a}], q[{b}];')
            elif e[0] == 'CCX':
                _, a, b, c = e
                lines.append(f'ccx q[{a}], q[{b}], q[{c}];')
            elif e[0] == 'H':
                _, a = e
                lines.append(f'h q[{a}];')
            elif e[0] == 'M':
                _, a = e
                lines.append(f'measure q[{a}] -> c[{m}];')
                m += 1
            elif e[0] == 'CPHASE':
                _, exponent, a, b = e
                if abs(exponent) >= cphase_angle_cutoff:
                    lines.append(f'cp(pi * {exponent}) q[{a}], q[{b}];')
            else:
                raise NotImplementedError(f'{e=}')
        return '\n'.join(lines)
