#include "flexible.h"
#include "gtest/gtest.h"

size_t N = 1000;
size_t M = 1000000;

template <typename Table>
void CommonTest() {
	MetaTable meta;

	const TypeID intType = meta.AddType("int", sizeof(int), __alignof(int));
	const TypeID doubleType = meta.AddType("double", sizeof(double), __alignof(double));
	const TypeID boolType = meta.AddType("bool", sizeof(bool), __alignof(bool));
	EXPECT_EQ(3, meta.GetTypeCount());

	const int iDefault = 1;
	const double dDefault = 2.0;
	const bool bDefault = true;
	const AttributeID iAttribute = meta.AddAttribute("i", intType, iDefault);
	const AttributeID dAttribute = meta.AddAttribute("d", doubleType, dDefault); 
	const AttributeID bAttribute = meta.AddAttribute("b", boolType, bDefault);
	EXPECT_EQ(3, meta.GetAttributeCount());

	Table t(meta);
	const Table& tc = t;	// constant version for testing const member functions

	EXPECT_EQ(&meta, &t.GetMetaTable());

	// Initial ReserveRows()
	t.ReserveRows(10);
	EXPECT_EQ(0, tc.GetRowCount());
	EXPECT_EQ(10, tc.GetRowCapacity());

	// Add new rows, which contains default values.
	t.AppendRows(5);
	EXPECT_EQ(5, tc.GetRowCount());
	for (size_t row = 0; row < tc.GetRowCount(); row++) {
		EXPECT_EQ(iDefault, tc.GetValue<int>(row, iAttribute));
		EXPECT_EQ(dDefault, tc.GetValue<double>(row, dAttribute));
		EXPECT_EQ(bDefault, tc.GetValue<bool>(row, bAttribute));
	}

	// Modify rows with SetValue()
	for (size_t row = 0; row < tc.GetRowCount(); row++) {
		int i = (int)row + 2;
		double d = row * 0.5;
		bool b = (row % 2) == 1;
		t.SetValue(row, iAttribute, i);
		t.SetValue(row, dAttribute, d);
		t.SetValue(row, bAttribute, b);
	}

	for (size_t row = 0; row < tc.GetRowCount(); row++) {
		int i = (int)row + 2;
		double d = row * 0.5;
		bool b = (row % 2) == 1;
		EXPECT_EQ(i, tc.GetValue<int>(row, iAttribute));
		EXPECT_EQ(d, tc.GetValue<double>(row, dAttribute));
		EXPECT_EQ(b, tc.GetValue<bool>(row, bAttribute));
	}

	// Copy constructor
	{
		Table t2(tc);
		EXPECT_EQ(&tc.GetMetaTable(), &t2.GetMetaTable());
		EXPECT_EQ(tc.GetRowCapacity(), t2.GetRowCapacity());
		EXPECT_EQ(tc.GetRowCount(), t2.GetRowCount());
		EXPECT_EQ(tc.GetMemorySize(), t2.GetMemorySize());

		for (size_t row = 0; row < tc.GetRowCount(); row++) {
			EXPECT_EQ(tc.GetValue<int>(row, iAttribute), t2.GetValue<int>(row, iAttribute));
			EXPECT_EQ(tc.GetValue<double>(row, dAttribute), t2.GetValue<double>(row, dAttribute));
			EXPECT_EQ(tc.GetValue<bool>(row, bAttribute), t2.GetValue<bool>(row, bAttribute));
		}

		// Change t2 (should not affect t as verified in later tests)
		for (size_t row = 0; row < t2.GetRowCount(); row++) {
			t2.SetValue(row, iAttribute, 0);
			t2.SetValue(row, dAttribute, 0.0);
			t2.SetValue(row, bAttribute, false);
		}
	}

	// Insert 3 rows()
	t.InsertRows(2, 3);
	EXPECT_EQ(5 + 3, tc.GetRowCount());
	for (size_t row = 0; row < tc.GetRowCount(); row++) {
		int i;
		double d;
		bool b;
		if (row < 2 || row >= 2 + 3) {
			int arow = row < 2 ? (int)row : (int)row - 3;
			i = arow + 2;
			d = arow * 0.5;
			b = (arow % 2) == 1;
		}
		else {
			i = iDefault;
			d = dDefault;
			b = bDefault;
		}

		EXPECT_EQ(i, tc.GetValue<int>(row, iAttribute));
		EXPECT_EQ(d, tc.GetValue<double>(row, dAttribute));
		EXPECT_EQ(b, tc.GetValue<bool>(row, bAttribute));
	}

	// Restore by RemoveRows()
	t.RemoveRows(2, 3);
	EXPECT_EQ(5, tc.GetRowCount());
	for (size_t row = 0; row < tc.GetRowCount(); row++) {
		int i = (int)row + 2;
		double d = row * 0.5;
		bool b = (row % 2) == 1;
		EXPECT_EQ(i, tc.GetValue<int>(row, iAttribute));
		EXPECT_EQ(d, tc.GetValue<double>(row, dAttribute));
		EXPECT_EQ(b, tc.GetValue<bool>(row, bAttribute));
	}

	// Reserve more to test reallocation
	t.ReserveRows(100);
	EXPECT_EQ(5, tc.GetRowCount());
	EXPECT_EQ(100, tc.GetRowCapacity());
	for (size_t row = 0; row < tc.GetRowCount(); row++) {
		int i = (int)row + 2;
		double d = row * 0.5;
		bool b = (row % 2) == 1;
		EXPECT_EQ(i, tc.GetValue<int>(row, iAttribute));
		EXPECT_EQ(d, tc.GetValue<double>(row, dAttribute));
		EXPECT_EQ(b, tc.GetValue<bool>(row, bAttribute));
	}

	std::cout << "Memory: " << tc.GetMemorySize() << std::endl;

	// Reserve less should make no changes
	t.ReserveRows(10);
	EXPECT_EQ(5, tc.GetRowCount());
	EXPECT_EQ(100, tc.GetRowCapacity());
}

TEST(AOSTable, Common) {
	CommonTest<AOSTable>();
}

TEST(SOATable, Common) {
	CommonTest<SOATable>();
}

TEST(FlexibleSOATable, Common) {
	CommonTest<FlexibleSOATable>();
}

void EulerIntegrationAOS(__m128* position, __m128* velocity, size_t stride) {
	const __m128 dt = _mm_set1_ps(0.01f);
	const __m128 a = _mm_setr_ps(0.0f, -9.8f, 0.0f, 0.0f);
	const __m128 adt = _mm_mul_ps(a, dt);

	for (size_t j = 0; j < M; j++) {
		__m128* __restrict p = position;
		__m128* __restrict v = velocity;

		for (size_t i = 0; i < N; i++) {
			*p = _mm_add_ps(*p, _mm_mul_ps(*v, dt));
			*v = _mm_add_ps(*v, adt);
			p = reinterpret_cast<__m128*>(reinterpret_cast<char*>(p) + stride);
			v = reinterpret_cast<__m128*>(reinterpret_cast<char*>(v) + stride);
		}
	}
}

TEST(AOS, EulerIntegration) {
	struct Particle {
		__m128 position;
		__m128 velocity;
	};

	Particle* particles = static_cast<Particle*>(_aligned_malloc(N * sizeof(Particle), __alignof(Particle)));

	for (size_t i = 0; i < N; i++) {
		particles[i].position = _mm_setr_ps(i + 0.0f, i + 1.0f, i + 2.0f, 0.0f);
		particles[i].velocity = _mm_setr_ps(i * 0.5f, i * 0.5f + 1.0f, i * 0.5f + 2.0f, 0.0f);
	}

	EulerIntegrationAOS(&particles->position, &particles->velocity, sizeof(Particle));

	// Check sum
	{
		float sum = 0.0f;
		for (Particle* p = particles; p != particles + N; p++) {
			float position[4];
			_mm_storeu_ps(position, p->position);
			sum += position[0] + position[1] + position[2];
		}
		std::cout << "CheckSum: " << sum << std::endl;
	}

	_aligned_free(particles);
}

TEST(AOSTable, EulerIntegration) {
	MetaTable meta;
	const TypeID vectorType = meta.AddType("vector", 16, 16);
	const AttributeID positionAttribute = meta.AddAttribute("position", vectorType, _mm_setzero_ps());
	const AttributeID velocityAttribute = meta.AddAttribute("velocity", vectorType, _mm_setzero_ps());

	AOSTable particles(meta);
	particles.ReserveRows(N);
	particles.AppendRows(N);
	std::cout << "Memory: " << particles.GetMemorySize() << std::endl;

	for (size_t i = 0; i < N; i++) {
		particles.SetValue(i, positionAttribute, _mm_setr_ps(i + 0.0f, i + 1.0f, i + 2.0f, 0.0f));
		particles.SetValue(i, velocityAttribute, _mm_setr_ps(i * 0.5f, i * 0.5f + 1.0f, i * 0.5f + 2.0f, 0.0f));
	}

	__m128* position = static_cast<__m128*>(particles.GetValueRaw(0, positionAttribute));
	__m128* velocity = static_cast<__m128*>(particles.GetValueRaw(0, velocityAttribute));	
	EulerIntegrationAOS(position, velocity, meta.GetAOSSize());

	// Check sum
	{
		float sum = 0.0f;
		const float* __restrict position = static_cast<const float*>(particles.GetValueRaw(0, positionAttribute));
		const size_t stride = meta.GetAOSSize();
		for (size_t i = 0; i < N; i++) {
			sum += position[0] + position[1] + position[2];
			position = reinterpret_cast<const float*>(reinterpret_cast<const char*>(position) + stride);
		}
		std::cout << "CheckSum: " << sum << std::endl;
	}
}

inline void StoreFloat3(float* p, const __m128 x) {
	_mm_store_ss(p, x);
	_mm_store_ss(p + 1, _mm_shuffle_ps(x, x, _MM_SHUFFLE(1, 1, 1, 1)));
	_mm_store_ss(p + 2, _mm_shuffle_ps(x, x, _MM_SHUFFLE(2, 2, 2, 2)));
}

void EulerIntegrationAOSFloat3(float* position, float* velocity, size_t stride) {
	const __m128 dt = _mm_set1_ps(0.01f);
	const __m128 a = _mm_setr_ps(0.0f, -9.8f, 0.0f, 0.0f);
	const __m128 adt = _mm_mul_ps(a, dt);

	for (size_t j = 0; j < M; j++) {
		float* __restrict p = position;
		float* __restrict v = velocity;

		for (size_t i = 0; i < N; i++) {
			const __m128 pi = _mm_loadu_ps(p);
			const __m128 vi = _mm_loadu_ps(v);
			StoreFloat3(p, _mm_add_ps(pi, _mm_mul_ps(vi, dt)));
			StoreFloat3(v, _mm_add_ps(vi, adt));
			p = reinterpret_cast<float*>(reinterpret_cast<char*>(p) + stride);
			v = reinterpret_cast<float*>(reinterpret_cast<char*>(v) + stride);
		}
	}
}

TEST(AOS, EulerIntegration_Float3) {
	struct Particle {
		float position[3];
		float velocity[3];
	};

	Particle* particles = static_cast<Particle*>(_aligned_malloc(N * sizeof(Particle), __alignof(Particle)));

	for (size_t i = 0; i < N; i++) {
		particles[i].position[0] = i + 0.0f;
		particles[i].position[1] = i + 1.0f;
		particles[i].position[2] = i + 2.0f;
		particles[i].velocity[0] = i * 0.5f + 0.0f;
		particles[i].velocity[1] = i * 0.5f + 1.0f;
		particles[i].velocity[2] = i * 0.5f + 2.0f;
	}

	EulerIntegrationAOSFloat3(particles->position, particles->velocity, sizeof(Particle));

	// Slower:
	//const __m128 dt = _mm_set1_ps(0.01f);
	//const __m128 a = _mm_setr_ps(0.0f, -9.8f, 0.0f, 0.0f);
	//const __m128 adt = _mm_mul_ps(a, dt);

	//for (size_t j = 0; j < M; j++) {
	//	for (Particle* p = particles; p != particles + N; p++) {
	//		const __m128 pi = _mm_loadu_ps(p->position);
	//		const __m128 vi = _mm_loadu_ps(p->velocity);
	//		StoreFloat3(p->position, _mm_add_ps(pi, _mm_mul_ps(vi, dt)));
	//		StoreFloat3(p->velocity, _mm_add_ps(vi, adt));
	//	}
	//}

	// Check sum
	{
		float sum = 0.0f;
		for (Particle* p = particles; p != particles + N; p++)
			sum += p->position[0] + p->position[1] + p->position[2];
		std::cout << "CheckSum: " << sum << std::endl;
	}

	_aligned_free(particles);
}

void EulerIntegrationAOSFloat3Naive(float* position, float* velocity, size_t stride) {
	const float dt = 0.01f;
	const float a[3] = { 0.0f, -9.8f, 0.0f };
	const float adt[3] = { a[0] * dt, a[1] * dt, a[2] * dt };

	for (size_t j = 0; j < M; j++) {
		float* __restrict p = position;
		float* __restrict v = velocity;

		for (size_t i = 0; i < N; i++) {
			p[0] += v[0] * dt;
			p[1] += v[1] * dt;
			p[2] += v[2] * dt;
			v[0] += adt[0];
			v[1] += adt[1];
			v[2] += adt[2];
			p = reinterpret_cast<float*>(reinterpret_cast<char*>(p) + stride);
			v = reinterpret_cast<float*>(reinterpret_cast<char*>(v) + stride);
		}
	}
}

TEST(AOS, EulerIntegration_Float3Naive) {
	struct Particle {
		float position[3];
		float velocity[3];
	};

	Particle* particles = static_cast<Particle*>(_aligned_malloc(N * sizeof(Particle), __alignof(Particle)));

	for (size_t i = 0; i < N; i++) {
		particles[i].position[0] = i + 0.0f;
		particles[i].position[1] = i + 1.0f;
		particles[i].position[2] = i + 2.0f;
		particles[i].velocity[0] = i * 0.5f + 0.0f;
		particles[i].velocity[1] = i * 0.5f + 1.0f;
		particles[i].velocity[2] = i * 0.5f + 2.0f;
	}

	EulerIntegrationAOSFloat3Naive(particles->position, particles->velocity, sizeof(Particle));

	// Check sum
	{
		float sum = 0.0f;
		for (Particle* p = particles; p != particles + N; p++)
			sum += p->position[0] + p->position[1] + p->position[2];
		std::cout << "CheckSum: " << sum << std::endl;
	}

	_aligned_free(particles);
}

TEST(AOSTable, EulerIntegration_Float3) {
	MetaTable meta;
	const TypeID float3Type = meta.AddType("float3", 12, 4);
	const float zero[3] = { 0.0f, 0.0f, 0.0f } ;
	const AttributeID positionAttribute = meta.AddAttribute("position", float3Type, zero);
	const AttributeID velocityAttribute = meta.AddAttribute("velocity", float3Type, zero);

	AOSTable particles(meta);
	particles.ReserveRows(N);
	particles.AppendRows(N);
	std::cout << "Memory: " << particles.GetMemorySize() << std::endl;

	for (size_t i = 0; i < N; i++) {
		const float p[3] = { i + 0.0f, i + 1.0f, i + 2.0f };
		const float v[3] = { i * 0.5f, i * 0.5f + 1.0f, i * 0.5f + 2.0f };
		particles.SetValue(i, positionAttribute, p);
		particles.SetValue(i, velocityAttribute, v);
	}

	const size_t stride = meta.GetAOSSize();
	float* __restrict position = static_cast<float*>(particles.GetValueRaw(0, positionAttribute));
	float* __restrict velocity = static_cast<float*>(particles.GetValueRaw(0, velocityAttribute));

	EulerIntegrationAOSFloat3(position, velocity, stride);

	// Check sum
	{
		float sum = 0.0f;
		const float* __restrict position = static_cast<const float*>(particles.GetValueRaw(0, positionAttribute));
		for (size_t i = 0; i < N; i++) {
			sum += position[0] + position[1] + position[2];
			position = reinterpret_cast<const float*>(reinterpret_cast<const char*>(position) + stride);
		}
		std::cout << "CheckSum: " << sum << std::endl;
	}
}

static void EulerIntegrationSOA(
	__m128* __restrict px,
	__m128* __restrict py,
	__m128* __restrict pz,
	__m128* __restrict vx,
	__m128* __restrict vy,
	__m128* __restrict vz)
{
	const __m128 dt = _mm_set1_ps(0.01f);
	const __m128 ax = _mm_set1_ps(0.0f);
	const __m128 ay = _mm_set1_ps(-9.8f);
	const __m128 az = _mm_set1_ps(0.0f);
	const __m128 axdt = _mm_mul_ps(ax, dt);
	const __m128 aydt = _mm_mul_ps(ay, dt);
	const __m128 azdt = _mm_mul_ps(az, dt);

	for (size_t j = 0; j < M; j++) {
		for (size_t i = 0; i < N / 4; i++) {
			px[i] = _mm_add_ps(px[i], _mm_mul_ps(vx[i], dt));
			py[i] = _mm_add_ps(py[i], _mm_mul_ps(vy[i], dt));
			pz[i] = _mm_add_ps(pz[i], _mm_mul_ps(vz[i], dt));
			vx[i] = _mm_add_ps(vx[i], axdt);
			vy[i] = _mm_add_ps(vy[i], aydt);
			vz[i] = _mm_add_ps(vz[i], azdt);
		}
	}
}

TEST(SOA, EulerIntegration) {
	float* positionX = static_cast<float*>(_aligned_malloc(N * 6 * sizeof(float), 16));
	float* positionY = positionX + N;
	float* positionZ = positionX + N * 2;
	float* velocityX = positionX + N * 3;
	float* velocityY = positionX + N * 4;
	float* velocityZ = positionX + N * 5;

	for (size_t i = 0; i < N; i++) {
		positionX[i] = i + 0.0f;
		positionY[i] = i + 1.0f;
		positionZ[i] = i + 2.0f;
		velocityX[i] = i * 0.5f;
		velocityY[i] = i * 0.5f + 1.0f;
		velocityZ[i] = i * 0.5f + 2.0f;
	}

	__m128* __restrict px = reinterpret_cast<__m128*>(positionX);
	__m128* __restrict py = reinterpret_cast<__m128*>(positionY);
	__m128* __restrict pz = reinterpret_cast<__m128*>(positionZ);
	__m128* __restrict vx = reinterpret_cast<__m128*>(velocityX);
	__m128* __restrict vy = reinterpret_cast<__m128*>(velocityY);
	__m128* __restrict vz = reinterpret_cast<__m128*>(velocityZ);

	EulerIntegrationSOA(px, py, pz, vx, vy, vz);

	// Check sum
	{
		float sum = 0.0f;

		for (size_t i = 0; i < N; i++)
			sum += positionX[i] + positionY[i] + positionZ[i];

		std::cout << "CheckSum: " << sum << std::endl;
	}

	_aligned_free(positionX);
}

TEST(SOATable, EulerIntegration) {
	MetaTable meta;
	const TypeID floatType = meta.AddType("float", 4, 16);
	const AttributeID positionXAttribute = meta.AddAttribute("positionX", floatType, 0.0f);
	const AttributeID positionYAttribute = meta.AddAttribute("positionY", floatType, 0.0f);
	const AttributeID positionZAttribute = meta.AddAttribute("positionZ", floatType, 0.0f);
	const AttributeID velocityXAttribute = meta.AddAttribute("velocityX", floatType, 0.0f);
	const AttributeID velocityYAttribute = meta.AddAttribute("velocityY", floatType, 0.0f);
	const AttributeID velocityZAttribute = meta.AddAttribute("velocityZ", floatType, 0.0f);

	SOATable particles(meta);
	particles.ReserveRows(N);
	particles.AppendRows(N);
	std::cout << "Memory: " << particles.GetMemorySize() << std::endl;

	for (size_t i = 0; i < N; i++) {
		particles.SetValue(i, positionXAttribute, i + 0.0f);
		particles.SetValue(i, positionYAttribute, i + 1.0f);
		particles.SetValue(i, positionZAttribute, i + 2.0f);
		particles.SetValue(i, velocityXAttribute, i * 0.5f);
		particles.SetValue(i, velocityYAttribute, i * 0.5f + 1.0f);
		particles.SetValue(i, velocityZAttribute, i * 0.5f + 2.0f);
	}

	__m128* __restrict px = static_cast<__m128*>(particles.GetValueRaw(0, positionXAttribute));
	__m128* __restrict py = static_cast<__m128*>(particles.GetValueRaw(0, positionYAttribute));
	__m128* __restrict pz = static_cast<__m128*>(particles.GetValueRaw(0, positionZAttribute));
	__m128* __restrict vx = static_cast<__m128*>(particles.GetValueRaw(0, velocityXAttribute));
	__m128* __restrict vy = static_cast<__m128*>(particles.GetValueRaw(0, velocityYAttribute));
	__m128* __restrict vz = static_cast<__m128*>(particles.GetValueRaw(0, velocityZAttribute));

	EulerIntegrationSOA(px, py, pz, vx, vy, vz);

	// Check sum
	{
		float sum = 0.0f;
		const float* __restrict px = static_cast<const float*>(particles.GetValueRaw(0, positionXAttribute));
		const float* __restrict py = static_cast<const float*>(particles.GetValueRaw(0, positionYAttribute));
		const float* __restrict pz = static_cast<const float*>(particles.GetValueRaw(0, positionZAttribute));
		for (size_t i = 0; i < N; i++)
			sum += px[i] + py[i] + pz[i];

		std::cout << "CheckSum: " << sum << std::endl;
	}
}

TEST(AOSTable, MinimumDistance_Naive) {
	MetaTable meta;
	const TypeID vectorType = meta.AddType("vector", 16, 16);
	const AttributeID sphereAttribute = meta.AddAttribute("sphere", vectorType, _mm_setzero_ps());	// XYZR

	AOSTable particles(meta);
	particles.ReserveRows(N);
	particles.AppendRows(N);
	std::cout << "Memory: " << particles.GetMemorySize() << std::endl;

	for (size_t i = 0; i < N; i++)
		particles.SetValue(i, sphereAttribute, _mm_setr_ps(i + 0.0f, i + 1.0f, i + 2.0f, i * 0.1f));

	const size_t stride = meta.GetAOSSize();
	float sum = 0.0f;
	for (size_t j = 0; j < M; j++) {
		float* __restrict sphere = static_cast<float*>(particles.GetValueRaw(0, sphereAttribute));
		const float q[3] = { j + 1.0f, 2.0f, 3.0f };
		float minD = FLT_MAX;

		for (size_t i = 0; i < N; i++) {
			const float *s = sphere;
			const float d[3] = { s[0] - q[0], s[1] - q[1], s[2] - q[2] };
			const float dd[3] = { d[0] * d[0], d[1] * d[1], d[2] * d[2] };
			const float dotxyz = dd[0] + dd[1] + dd[2];
			const float r = s[3];
			const float sd = sqrtf(dotxyz) - r;
			minD = std::min(minD, sd);
			sphere = reinterpret_cast<float*>(reinterpret_cast<char*>(sphere) + stride);
		}
		sum += minD;
	}

	std::cout << "CheckSum: " << sum << std::endl;
}

TEST(AOSTable, MinimumDistance) {
	MetaTable meta;
	const TypeID vectorType = meta.AddType("vector", 16, 16);
	const AttributeID sphereAttribute = meta.AddAttribute("sphere", vectorType, _mm_setzero_ps());	// XYZR

	AOSTable particles(meta);
	particles.ReserveRows(N);
	particles.AppendRows(N);
	std::cout << "Memory: " << particles.GetMemorySize() << std::endl;

	for (size_t i = 0; i < N; i++)
		particles.SetValue(i, sphereAttribute, _mm_setr_ps(i + 0.0f, i + 1.0f, i + 2.0f, i * 0.1f));

	const size_t stride = meta.GetAOSSize();
	__m128 sum = _mm_setzero_ps();
	for (size_t j = 0; j < M; j++) {
		__m128* __restrict sphere = static_cast<__m128*>(particles.GetValueRaw(0, sphereAttribute));
		const __m128 q = _mm_setr_ps(j + 1.0f, 2.0f, 3.0f, 0.0f);
		__m128 minD = _mm_set1_ps(FLT_MAX);

		for (size_t i = 0; i < N; i++) {
			const __m128 s = *sphere;
			const __m128 d = _mm_sub_ps(s, q);
			const __m128 dd = _mm_mul_ps(d, d);
			const __m128 dotxy = _mm_add_ss(dd, _mm_shuffle_ps(dd, dd, _MM_SHUFFLE(1, 1, 1, 1)));
			const __m128 dotxyz = _mm_add_ss(dotxy, _mm_shuffle_ps(dd, dd, _MM_SHUFFLE(2, 2, 2, 2)));
			const __m128 r = _mm_shuffle_ps(s, s, _MM_SHUFFLE(3, 3, 3, 3));
			const __m128 sd = _mm_sub_ss(_mm_sqrt_ss(dotxyz), r);
			minD = _mm_min_ss(minD, sd);
			sphere = reinterpret_cast<__m128*>(reinterpret_cast<char*>(sphere) + stride);
		}
		sum = _mm_add_ss(sum, minD);
	}

	std::cout << "CheckSum: " << _mm_cvtss_f32(sum) << std::endl;
}

TEST(SOATable, MinimumDistance) {
	MetaTable meta;
	const TypeID vectorType = meta.AddType("float", 4, 16);
	const AttributeID sphereXAttribute = meta.AddAttribute("sphereX", vectorType, 0.0f);
	const AttributeID sphereYAttribute = meta.AddAttribute("sphereY", vectorType, 0.0f);
	const AttributeID sphereZAttribute = meta.AddAttribute("sphereZ", vectorType, 0.0f);
	const AttributeID sphereRAttribute = meta.AddAttribute("sphereR", vectorType, 0.0f);

	SOATable particles(meta);
	particles.ReserveRows(N);
	particles.AppendRows(N);
	std::cout << "Memory: " << particles.GetMemorySize() << std::endl;

	for (size_t i = 0; i < N; i++) {
		particles.SetValue(i, sphereXAttribute, i + 0.0f);
		particles.SetValue(i, sphereYAttribute, i + 1.0f);
		particles.SetValue(i, sphereZAttribute, i + 2.0f);
		particles.SetValue(i, sphereRAttribute, i * 0.1f);
	}

	const size_t stride = meta.GetAOSSize();
	__m128 sum = _mm_setzero_ps();
	for (size_t j = 0; j < M; j++) {
		__m128* __restrict sx = static_cast<__m128*>(particles.GetValueRaw(0, sphereXAttribute));
		__m128* __restrict sy = static_cast<__m128*>(particles.GetValueRaw(0, sphereYAttribute));
		__m128* __restrict sz = static_cast<__m128*>(particles.GetValueRaw(0, sphereZAttribute));
		__m128* __restrict sr = static_cast<__m128*>(particles.GetValueRaw(0, sphereRAttribute));
		const __m128 qx = _mm_set1_ps(j + 1.0f);
		const __m128 qy = _mm_set1_ps(2.0f);
		const __m128 qz = _mm_set1_ps(3.0f);
		__m128 minD = _mm_set1_ps(FLT_MAX);

		for (size_t i = 0; i < N / 4; i++) {
			const __m128 dx = _mm_sub_ps(sx[i], qx);
			const __m128 dy = _mm_sub_ps(sy[i], qy);
			const __m128 dz = _mm_sub_ps(sz[i], qz);
			const __m128 dot = _mm_add_ps(_mm_add_ps(_mm_mul_ps(dx, dx), _mm_mul_ps(dy, dy)), _mm_mul_ps(dz, dz));
			const __m128 sd = _mm_sub_ps(_mm_sqrt_ps(dot), sr[i]);
			minD = _mm_min_ps(minD, sd);
		}
		minD = _mm_min_ps(minD, _mm_shuffle_ps(minD, minD, _MM_SHUFFLE(1, 0, 3, 2)));
		minD = _mm_min_ss(minD, _mm_shuffle_ps(minD, minD, _MM_SHUFFLE(1, 1, 1, 1)));

		sum = _mm_add_ss(sum, minD);
	}

	std::cout << "CheckSum: " << _mm_cvtss_f32(sum) << std::endl;
}

template<typename Table>
void OptimalSizeTest() {
	MetaTable meta;

	const TypeID a1s1t = meta.AddType("a1s1t", 1, 1);
	const TypeID a2s2t = meta.AddType("a2s2t", 2, 2);
	const TypeID a4s4t = meta.AddType("a4s4t", 4, 4);

	meta.AddAttribute("a1s1", a1s1t, (char)0);
	meta.AddAttribute("a2s2", a2s2t, (short)0);
	meta.AddAttribute("a4s4", a4s4t, (int)0);
	
	Table table(meta);

	table.ReserveRows(1);

	EXPECT_EQ(8, table.GetAttributeSize());
}

TEST(FlexibleSOATable, OptimalSize) {
	OptimalSizeTest<FlexibleSOATable>();
}

TEST(SOATable, OptimalSize) {
	OptimalSizeTest<SOATable>();
}

TEST(FlexibleSOATable, Variability) {
	MetaTable meta;

	const TypeID intType = meta.AddType("int", sizeof(int), __alignof(int));
	const TypeID doubleType = meta.AddType("double", sizeof(double), __alignof(double));
	const TypeID boolType = meta.AddType("bool", sizeof(bool), __alignof(bool));

	const int iDefault = 1;
	const double dDefault = 2.0;
	const bool bDefault = true;
	const AttributeID iAttribute = meta.AddAttribute("i", intType, iDefault);
	const AttributeID dAttribute = meta.AddAttribute("d", doubleType, dDefault); 
	const AttributeID bAttribute = meta.AddAttribute("b", boolType, bDefault);

	FlexibleSOATable t(meta);

	// Even empty table can contains uniform values.
	for (AttributeID attribute = 0; attribute != meta.GetAttributeCount(); attribute++)
		EXPECT_TRUE(t.IsColumnUniform(attribute));

	// Add new rows, which contains default values.
	t.AppendRows(5);
	EXPECT_EQ(5, t.GetRowCount());

	// Still be uniform
	for (AttributeID attribute = 0; attribute != meta.GetAttributeCount(); attribute++)
		EXPECT_TRUE(t.IsColumnUniform(attribute));

	// Use GetUniformValue()
	EXPECT_EQ(iDefault, t.GetUniformValue<int>(iAttribute));
	EXPECT_EQ(dDefault, t.GetUniformValue<double>(dAttribute));
	EXPECT_EQ(bDefault, t.GetUniformValue<bool>(bAttribute));

	// Use GetValue() with const
	const FlexibleSOATable& tc = t;
	for (size_t row = 0; row < t.GetRowCount(); row++) {
		EXPECT_EQ(iDefault, tc.GetValue<int>(row, iAttribute));
		EXPECT_EQ(dDefault, tc.GetValue<double>(row, dAttribute));
		EXPECT_EQ(bDefault, tc.GetValue<bool>(row, bAttribute));
	}
	
	// Still be uniform
	for (AttributeID attribute = 0; attribute != meta.GetAttributeCount(); attribute++)
		EXPECT_TRUE(t.IsColumnUniform(attribute));

	// Change i to varying by SetValue()
	for (size_t row = 0; row < t.GetRowCount(); row++)
		t.SetValue(row, iAttribute, (int)row + 1);

	// Changed to varying but not affecting other columns
	EXPECT_FALSE(t.IsColumnUniform(iAttribute));
	EXPECT_TRUE(t.IsColumnUniform(dAttribute));
	EXPECT_TRUE(t.IsColumnUniform(bAttribute));

	for (size_t row = 0; row < t.GetRowCount(); row++) {
		EXPECT_EQ(row + 1, tc.GetValue<int>(row, iAttribute));
		EXPECT_EQ(row + 1, t.GetValue<int>(row, iAttribute));
	}

	// GetVaryingValues()
	{
		double* d = static_cast<double*>(t.GetVaryingValuesRaw(dAttribute));
		
		EXPECT_FALSE(t.IsColumnUniform(dAttribute));

		for (size_t row = 0; row < t.GetRowCount(); row++)
			EXPECT_EQ(dDefault, d[row]);

		// Modify values in array directly
		for (size_t row = 0; row < t.GetRowCount(); row++)
			d[row] = row + 1.0;
	}

	for (size_t row = 0; row < t.GetRowCount(); row++) {
		EXPECT_EQ(row + 1.0, tc.GetValue<double>(row, dAttribute));
		EXPECT_EQ(row + 1.0, t.GetValue<double>(row, dAttribute));
	}

	// Change back to uniform
	{
		double d = -1;
		t.SetUniformValue(dAttribute, d);

		EXPECT_TRUE(tc.IsColumnUniform(dAttribute));
		EXPECT_EQ(d, t.GetUniformValue<double>(dAttribute));
	}
}

int main(int argc, char* argv[]) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
