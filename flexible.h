#pragma once

#include <cassert>
#include <string>
#include <vector>

////////////////////////////////////////////////////////////////////////////////

typedef size_t TypeID;
typedef size_t AttributeID;

struct MetaType {
	std::string name;
	size_t size;
	size_t alignment;
};

struct MetaAttribute {
	std::string name;
	TypeID type;
	void* defaultValue;
	size_t AOSOffset;
};

class MetaTable {
public:
	MetaTable();
	~MetaTable();

	size_t GetTypeCount() const;
	const MetaType& GetType(TypeID type) const;
	TypeID FindType(const char* name) const;
	TypeID AddType(const char* name, size_t size, size_t alignment);

	size_t GetAttributeCount() const;
	const MetaAttribute& GetAttribute(AttributeID attribute) const;
	AttributeID FindAttribute(const char* name);

	AttributeID AddAttributeRaw(const char* name, TypeID type, const void* defaultValue);

	template<typename T> AttributeID AddAttribute(const char* name, TypeID type, const T& defaultValue);

	size_t GetAOSSize() const;
	size_t GetAOSAlignment() const;

private:
	typedef std::vector<MetaType> MetaTypeList;
	typedef std::vector<MetaAttribute> MetaAttributeList;

	MetaTypeList mTypes;
	MetaAttributeList mAttributes;

	size_t mAOSSize;
	size_t mAOSAlignment;

	bool mMutable;
};

inline MetaTable::MetaTable() : mAOSSize(0), mAOSAlignment(0) {
}

inline MetaTable::~MetaTable() {
	for (MetaAttributeList::iterator itr = mAttributes.begin(); itr != mAttributes.end(); ++itr)
		_aligned_free(itr->defaultValue);
}

inline size_t MetaTable::GetTypeCount() const {
	return mTypes.size();
}

inline const MetaType& MetaTable::GetType(TypeID type) const {
	return mTypes[type];
}

inline TypeID MetaTable::FindType(const char* name) const {
	assert(name != nullptr);

	for (size_t i = 0; i < mTypes.size(); i++)
		if (mTypes[i].name == name)
			return i;

	return static_cast<size_t>(-1);
}

inline TypeID MetaTable::AddType(const char* name, size_t size, size_t alignment) {
	assert(FindType(name) == static_cast<size_t>(-1));
	assert(size >= 1);
	assert(alignment != 0 && (alignment & (~alignment + 1)) == alignment);	// Alignment must be power of 2.

	MetaType t;
	t.name = name;
	t.size = size;
	t.alignment = alignment;

	TypeID typeID = mTypes.size();
	mTypes.push_back(t);
	return typeID;
}

inline size_t MetaTable::GetAttributeCount() const {
	return mAttributes.size();
}

inline const MetaAttribute& MetaTable::GetAttribute(AttributeID attribute) const {
	return mAttributes[attribute];
}

inline AttributeID MetaTable::FindAttribute(const char* name) {
	assert(name != nullptr);

	for (size_t i = 0; i < mAttributes.size(); i++)
		if (mAttributes[i].name == name)
			return i;

	return static_cast<size_t>(-1);
}

inline AttributeID MetaTable::AddAttributeRaw(const char* name, TypeID type, const void* defaultValue) {
	assert(FindAttribute(name) == static_cast<size_t>(-1));
	assert(defaultValue != nullptr);

	const MetaType& t = GetType(type);

	MetaAttribute a;
	a.name = name;
	a.type = type;
	a.defaultValue = _aligned_malloc(t.size, t.alignment);
	memcpy(a.defaultValue, defaultValue, t.size);

	// Compute AOS information
	if (mAOSAlignment < t.alignment)
		mAOSAlignment = t.alignment;

	a.AOSOffset = mAttributes.empty() ? 0 : (mAttributes.back().AOSOffset + GetType(mAttributes.back().type).size + t.alignment - 1) & ~(t.alignment - 1);
	mAOSSize = (a.AOSOffset + t.size + mAOSAlignment - 1) & ~(mAOSAlignment - 1);

	const AttributeID attributeID = mAttributes.size();
	mAttributes.push_back(a);
	return attributeID;
}

template<typename T>
AttributeID MetaTable::AddAttribute(const char* name, TypeID type, const T& defaultValue) {
	assert(sizeof(T) == GetType(type).size);
	return AddAttributeRaw(name, type, &defaultValue);
}

inline size_t MetaTable::GetAOSSize() const {
	return mAOSSize;
}

inline size_t MetaTable::GetAOSAlignment() const {
	return mAOSAlignment;
}

////////////////////////////////////////////////////////////////////////////////

class AOSTable {
public:
	explicit AOSTable(const MetaTable& meta);
	AOSTable(const AOSTable& rhs);
	~AOSTable();

	const MetaTable& GetMetaTable() const;
	
	size_t GetRowCapacity() const;
	size_t GetRowCount() const;

	void ReserveRows(size_t rowCapacity);
	void AppendRows(size_t rowCount);
	void InsertRows(size_t rowIndex, size_t rowCount);
	void RemoveRows(size_t rowIndex, size_t rowCount);

	const void* GetValueRaw(size_t rowIndex, AttributeID attribute) const;
	void* GetValueRaw(size_t rowIndex, AttributeID attribute);
	void SetValueRaw(size_t rowIndex, AttributeID attribute, const void* value);

	template <typename T> const T& GetValue(size_t rowIndex, AttributeID attribute) const;
	template <typename T> T& GetValue(size_t rowIndex, AttributeID attribute);
	template <typename T> void SetValue(size_t rowIndex, AttributeID attribute, const T& value);

	size_t GetMemorySize() const;

private:
	void ReserveRowsInternal(size_t rowCount);

	const MetaTable& mMeta;
	void *mRowData;
	size_t mRowCapacity;
	size_t mRowCount;
};

inline AOSTable::AOSTable(const MetaTable& meta) : 
	mMeta(meta),
	mRowData(nullptr),
	mRowCapacity(0),
	mRowCount(0)
{
}

inline AOSTable::AOSTable(const AOSTable& rhs) : 
	mMeta(rhs.mMeta), 
	mRowData(nullptr),
	mRowCapacity(rhs.mRowCapacity),
	mRowCount(rhs.mRowCount)
{
	mRowData = _aligned_malloc(mMeta.GetAOSSize() * mRowCapacity, mMeta.GetAOSAlignment());
	memcpy(mRowData, rhs.mRowData, mMeta.GetAOSSize() * mRowCount);	
}

inline AOSTable::~AOSTable() {
	_aligned_free(mRowData);
}

inline const MetaTable& AOSTable::GetMetaTable() const {
	return mMeta;
}

inline size_t AOSTable::GetRowCapacity() const {
	return mRowCapacity;
}

inline size_t AOSTable::GetRowCount() const {
	return mRowCount;
}

inline void AOSTable::ReserveRows(size_t rowCapacity) {
	if (rowCapacity > mRowCapacity) {
		void *rowData = _aligned_malloc(mMeta.GetAOSSize() * rowCapacity, mMeta.GetAOSAlignment());
		memcpy(rowData, mRowData, mMeta.GetAOSSize()  * mRowCount);
		_aligned_free(mRowData);
		mRowData = rowData;
		mRowCapacity = rowCapacity;
	}
}

inline void AOSTable::ReserveRowsInternal(size_t rowCount) {
	if (mRowCapacity < rowCount) {
		size_t newCapacity = mRowCapacity == 0 ? 16 : mRowCapacity * 2;
		while (newCapacity < rowCount)
			newCapacity *= 2;
		ReserveRows(newCapacity);
	}
}

inline void AOSTable::AppendRows(size_t rowCount) {
	const size_t newCount = mRowCount + rowCount;
	ReserveRowsInternal(newCount);

	// Initialize new rows
	size_t rowIndex = mRowCount;
	mRowCount = newCount; // Need to do this before calling SetValue()
	for (; rowIndex < newCount; rowIndex++)
		for (AttributeID attribute = 0; attribute != mMeta.GetAttributeCount(); attribute++)
			SetValueRaw(rowIndex, attribute, mMeta.GetAttribute(attribute).defaultValue);
}

inline void AOSTable::InsertRows(size_t rowIndex, size_t rowCount) {
	assert(rowIndex < mRowCount);
	
	// To be optimized: reallocate and move first/second half.
	const size_t newCount = mRowCount + rowCount;
	ReserveRowsInternal(newCount);

	// Move second half to the back.
	const size_t size = mMeta.GetAOSSize();
	char* p = static_cast<char*>(mRowData) + rowIndex * size;
	memmove(p + rowCount * size, p, (mRowCount - rowIndex) * size);

	// Initialize new rows
	const size_t end = rowIndex + rowCount;
	mRowCount = newCount; // Need to do this before calling SetValue()
	for (; rowIndex < end; rowIndex++)
		for (AttributeID attribute = 0; attribute != mMeta.GetAttributeCount(); attribute++)
			SetValueRaw(rowIndex, attribute, mMeta.GetAttribute(attribute).defaultValue);
}

inline void AOSTable::RemoveRows(size_t rowIndex, size_t rowCount) {
	assert(rowIndex + rowCount <= mRowCount);
	
	if (rowCount > 0) {
		// move second half
		const size_t size = mMeta.GetAOSSize();
		char* p = static_cast<char*>(mRowData) + rowIndex * size;
		memmove(p, p + rowCount * size, (mRowCount - rowIndex - rowCount) * size);
	
		mRowCount -= rowCount;
	}
}

inline const void* AOSTable::GetValueRaw(size_t rowIndex, AttributeID attribute) const {
	return const_cast<AOSTable *>(this)->GetValueRaw(rowIndex, attribute);
}

inline void* AOSTable::GetValueRaw(size_t rowIndex, AttributeID attribute) {
	assert(rowIndex < mRowCount);
	const MetaAttribute& a = mMeta.GetAttribute(attribute);
	return static_cast<char*>(mRowData) + rowIndex * mMeta.GetAOSSize() + a.AOSOffset;
}

inline void AOSTable::SetValueRaw(size_t rowIndex, AttributeID attribute, const void* value) {
	assert(rowIndex < mRowCount);
	const MetaAttribute& a = mMeta.GetAttribute(attribute);
	const MetaType& t = mMeta.GetType(a.type);
	memcpy(static_cast<char*>(mRowData) + rowIndex * mMeta.GetAOSSize() + a.AOSOffset, value, t.size);
}

template <typename T>
const T& AOSTable::GetValue(size_t rowIndex, AttributeID attribute) const {
	assert(sizeof(T) == mMeta.GetType(mMeta.GetAttribute(attribute).type).size);
	return static_cast<const T*>(GetValueRaw(rowIndex, attribute))[0];
}

template <typename T>
T& AOSTable::GetValue(size_t rowIndex, AttributeID attribute) {
	assert(sizeof(T) == mMeta.GetType(mMeta.GetAttribute(attribute).type).size);
	return static_cast<T*>(GetValueRaw(rowIndex, attribute))[0];
}

template <typename T>
void AOSTable::SetValue(size_t rowIndex, AttributeID attribute, const T& value) {
	assert(sizeof(T) == mMeta.GetType(mMeta.GetAttribute(attribute).type).size);
	SetValueRaw(rowIndex, attribute, &value);
}

inline size_t AOSTable::GetMemorySize() const {
	return sizeof(this) 
		+ mRowCapacity * mMeta.GetAOSSize();		// mRowData
}

////////////////////////////////////////////////////////////////////////////////

class SOATable {
public:
	explicit SOATable(const MetaTable& meta);
	SOATable(const SOATable& rhs);
	~SOATable();

	const MetaTable& GetMetaTable() const;
	
	size_t GetRowCapacity() const;
	size_t GetRowCount() const;

	void ReserveRows(size_t rowCapacity);
	void AppendRows(size_t rowCount);
	void InsertRows(size_t rowIndex, size_t rowCount);
	void RemoveRows(size_t rowIndex, size_t rowCount);

	const void* GetValueRaw(size_t rowIndex, AttributeID attribute) const;
	void* GetValueRaw(size_t rowIndex, AttributeID attribute);
	void SetValueRaw(size_t rowIndex, AttributeID attribute, const void* value);

	template <typename T> const T& GetValue(size_t rowIndex, AttributeID attribute) const;
	template <typename T> T& GetValue(size_t rowIndex, AttributeID attribute);
	template <typename T> void SetValue(size_t rowIndex, AttributeID attribute, const T& value);

	size_t GetAttributeSize() const;
	size_t GetMemorySize() const;

private:
	void ReserveRowsInternal(size_t rowCount);

	const MetaTable& mMeta;
	void **mColumnData;
	size_t mRowCapacity;
	size_t mRowCount;
};

inline SOATable::SOATable(const MetaTable& meta) : 
	mMeta(meta),
	mColumnData(nullptr),
	mRowCapacity(0),
	mRowCount(0)
{
	assert(mMeta.GetAttributeCount() > 0);

	mColumnData = (void**)malloc(mMeta.GetAttributeCount() * sizeof(void*));
	for (AttributeID attribute = 0; attribute != mMeta.GetAttributeCount(); attribute++)
		mColumnData[attribute] = nullptr;
}

inline SOATable::SOATable(const SOATable& rhs) :
	mMeta(rhs.mMeta),
	mColumnData(nullptr),
	mRowCapacity(0),
	mRowCount(0)
{
	assert(mMeta.GetAttributeCount() > 0);

	ReserveRows(rhs.mRowCapacity);
	for (AttributeID attribute = 0; attribute != mMeta.GetAttributeCount(); attribute++) {
		const MetaType& t = mMeta.GetType(mMeta.GetAttribute(attribute).type);
		memcpy(mColumnData[attribute], rhs.mColumnData[attribute], t.size * rhs.mRowCount);
	}
	mRowCapacity = rhs.mRowCapacity;
	mRowCount = rhs.mRowCount;
}

inline SOATable::~SOATable() {
	_aligned_free(mColumnData[0]);
	free(mColumnData);
}

inline const MetaTable& SOATable::GetMetaTable() const {
	return mMeta;
}

inline size_t SOATable::GetRowCapacity() const {
	return mRowCapacity;
}

inline size_t SOATable::GetRowCount() const {
	return mRowCount;
}

inline void SOATable::ReserveRows(size_t rowCapacity) {
	if (rowCapacity > mRowCapacity) {
		size_t newSize = 0;
		void** columnData = (void**)malloc(mMeta.GetAttributeCount() * sizeof(void*));
		for (AttributeID attribute = 0; attribute != mMeta.GetAttributeCount(); attribute++) {
			const MetaType& t = mMeta.GetType(mMeta.GetAttribute(attribute).type);
			// newSize = (newSize + mMeta.GetAOSAlignment() - 1) & ~(mMeta.GetAOSAlignment() - 1);
			newSize = (newSize + t.alignment - 1) & ~(t.alignment - 1);
			columnData[attribute] = (void*)newSize;	// Store offset first
			newSize += t.size * rowCapacity;
		}

		size_t fixup = (size_t)_aligned_malloc(newSize, mMeta.GetAOSAlignment());

		for (AttributeID attribute = 0; attribute != mMeta.GetAttributeCount(); attribute++) {
			const MetaType& t = mMeta.GetType(mMeta.GetAttribute(attribute).type);
			columnData[attribute] = static_cast<char*>(columnData[attribute]) + fixup;
			if (mColumnData)
				memcpy(columnData[attribute], mColumnData[attribute], t.size * mRowCount);
		}

		if (mColumnData) {
			_aligned_free(mColumnData[0]);
			free(mColumnData);
		}
		
		mColumnData = columnData;
		mRowCapacity = rowCapacity;
	}
}

inline void SOATable::ReserveRowsInternal(size_t rowCount) {
	if (mRowCapacity < rowCount) {
		size_t newCapacity = mRowCapacity == 0 ? 16 : mRowCapacity * 2;
		while (newCapacity < rowCount)
			newCapacity *= 2;
		ReserveRows(newCapacity);
	}
}

inline void SOATable::AppendRows(size_t rowCount) {
	size_t newCount = mRowCount + rowCount;
	ReserveRowsInternal(rowCount);

	// Initialize new rows
	size_t start = mRowCount;
	mRowCount = newCount; // Need to do this before calling SetValue()
	for (AttributeID attribute = 0; attribute != mMeta.GetAttributeCount(); attribute++) {
		const void* defaultValue = mMeta.GetAttribute(attribute).defaultValue;
		for (size_t rowIndex = start; rowIndex < newCount; rowIndex++)
			SetValueRaw(rowIndex, attribute, defaultValue);
	}
}

inline void SOATable::InsertRows(size_t rowIndex, size_t rowCount) {
	assert(rowIndex < mRowCount);

	// To be optimized: reallocate and move first/second half.
	const size_t newCount = mRowCount + rowCount;
	ReserveRowsInternal(newCount);

	// Move second half to the back.
	for (AttributeID attribute = 0; attribute != mMeta.GetAttributeCount(); attribute++) {
		const size_t size = mMeta.GetType(mMeta.GetAttribute(attribute).type).size;
		char* p = static_cast<char*>(mColumnData[attribute]) + rowIndex * size;
		memmove(p + rowCount * size, p, (mRowCount - rowIndex) * size);
	}

	// Initialize new rows
	const size_t start = rowIndex;
	const size_t end = rowIndex + rowCount;
	mRowCount = newCount; // Need to do this before calling SetValue()
	for (AttributeID attribute = 0; attribute != mMeta.GetAttributeCount(); attribute++) {
		const void* defaultValue = mMeta.GetAttribute(attribute).defaultValue;
		for (size_t i = start; i != end; ++i)
			SetValueRaw(i, attribute, defaultValue);
	}
}

inline void SOATable::RemoveRows(size_t rowIndex, size_t rowCount) {
	assert(rowIndex + rowCount <= mRowCount);

	if (rowCount > 0) {
		// move second half
		for (AttributeID attribute = 0; attribute != mMeta.GetAttributeCount(); attribute++) {
			const size_t size = mMeta.GetType(mMeta.GetAttribute(attribute).type).size;
			char* p = static_cast<char*>(mColumnData[attribute]) + rowIndex * size;
			memmove(p, p + rowCount * size, (mRowCount - rowIndex - rowCount) * size);
		}

		mRowCount -= rowCount;
	}
}

inline const void* SOATable::GetValueRaw(size_t rowIndex, AttributeID attribute) const {
	return const_cast<SOATable *>(this)->GetValueRaw(rowIndex, attribute);
}

inline void* SOATable::GetValueRaw(size_t rowIndex, AttributeID attribute) {
	const MetaAttribute& a = mMeta.GetAttribute(attribute);
	const MetaType& t = mMeta.GetType(a.type);
	return static_cast<char*>(mColumnData[attribute]) + rowIndex * t.size;
}

inline void SOATable::SetValueRaw(size_t rowIndex, AttributeID attribute, const void* value) {
	const MetaAttribute& a = mMeta.GetAttribute(attribute);
	const MetaType& t = mMeta.GetType(a.type);
	memcpy(static_cast<char*>(mColumnData[attribute]) + rowIndex * t.size, value, t.size);
}

template <typename T>
const T& SOATable::GetValue(size_t rowIndex, AttributeID attribute) const {
	assert(sizeof(T) == mMeta.GetType(mMeta.GetAttribute(attribute).type).size);
	return static_cast<const T*>(GetValueRaw(rowIndex, attribute))[0];
}

template <typename T>
T& SOATable::GetValue(size_t rowIndex, AttributeID attribute) {
	assert(sizeof(T) == mMeta.GetType(mMeta.GetAttribute(attribute).type).size);
	return static_cast<T*>(GetValueRaw(rowIndex, attribute))[0];
}

template <typename T>
void SOATable::SetValue(size_t rowIndex, AttributeID attribute, const T& value) {
	assert(sizeof(T) == mMeta.GetType(mMeta.GetAttribute(attribute).type).size);
	SetValueRaw(rowIndex, attribute, &value);
}

inline size_t SOATable::GetAttributeSize() const {
	size_t newSize = 0;
	for (AttributeID attribute = 0; attribute != mMeta.GetAttributeCount(); attribute++) {
		const MetaType& t = mMeta.GetType(mMeta.GetAttribute(attribute).type);
		// newSize = (newSize + mMeta.GetAOSAlignment() - 1) & ~(mMeta.GetAOSAlignment() - 1);
		newSize = (newSize + t.alignment - 1) & ~(t.alignment - 1);
		newSize += t.size * mRowCapacity;
	}
	return newSize;
}

inline size_t SOATable::GetMemorySize() const {
	return GetAttributeSize()
		+ sizeof(this) 
		+ mMeta.GetAttributeCount() * sizeof(void*);// mColumnData
}

////////////////////////////////////////////////////////////////////////////////

// Added variability attributes with specialized API
class FlexibleSOATable {
public:
	explicit FlexibleSOATable(const MetaTable& meta);
	FlexibleSOATable(const FlexibleSOATable& rhs);
	~FlexibleSOATable();

	const MetaTable& GetMetaTable() const;
	
	size_t GetRowCapacity() const;
	size_t GetRowCount() const;

	void ReserveRows(size_t rowCapacity);
	void AppendRows(size_t rowCount);
	void InsertRows(size_t rowIndex, size_t rowCount);
	void RemoveRows(size_t rowIndex, size_t rowCount);

	const void* GetValueRaw(size_t rowIndex, AttributeID attribute) const;
	void* GetValueRaw(size_t rowIndex, AttributeID attribute);
	void SetValueRaw(size_t rowIndex, AttributeID attribute, const void* value);

	template <typename T> const T& GetValue(size_t rowIndex, AttributeID attribute) const;
	template <typename T> T& GetValue(size_t rowIndex, AttributeID attribute);
	template <typename T> void SetValue(size_t rowIndex, AttributeID attribute, const T& value);

	size_t GetAttributeSize() const;
	size_t GetMemorySize() const;

	// Variability attributes

	bool IsColumnUniform(AttributeID attribute) const;

	const void* GetUniformValueRaw(AttributeID attribute) const;
	void SetUniformValueRaw(AttributeID attribute, const void* value);

	template<typename T> const T& GetUniformValue(AttributeID attribute) const;
	template<typename T> void SetUniformValue(AttributeID attribute, const T& value);

	void* GetVaryingValuesRaw(AttributeID attribute);
	void SetVaryingValuesRaw(size_t rowIndex, size_t rowCount, AttributeID attribute, const void* values);

private:
	void ReserveRowsInternal(size_t rowCount);
	void ConvertColumnToVarying(AttributeID attribute);

	struct Column {
		void* mData;		// [0] is uniform, [0..mRowCapacity-1] is varying
		bool mUniform;
	};

	const MetaTable& mMeta;
	Column *mColumns;
	size_t mRowCapacity;
	size_t mRowCount;
};

inline FlexibleSOATable::FlexibleSOATable(const MetaTable& meta) : 
	mMeta(meta),
	mColumns(nullptr),
	mRowCapacity(1),	// For storing uniform
	mRowCount(0)
{
	assert(mMeta.GetAttributeCount() > 0);

	mColumns = (Column*)malloc(mMeta.GetAttributeCount() * sizeof(Column));

	size_t newSize = 0;
	for (AttributeID attribute = 0; attribute != mMeta.GetAttributeCount(); attribute++) {
		const MetaType& t = mMeta.GetType(mMeta.GetAttribute(attribute).type);
		newSize = (newSize + mMeta.GetAOSAlignment() - 1) & ~(mMeta.GetAOSAlignment() - 1);
		mColumns[attribute].mData = (void*)newSize;
		mColumns[attribute].mUniform = true;
		newSize += t.size * mRowCapacity;
	}

	size_t fixup = (size_t)_aligned_malloc(newSize, mMeta.GetAOSAlignment());

	for (AttributeID attribute = 0; attribute != mMeta.GetAttributeCount(); attribute++) {
		const MetaAttribute& a = mMeta.GetAttribute(attribute);
		const MetaType& t = mMeta.GetType(a.type);
		mColumns[attribute].mData = static_cast<char*>(mColumns[attribute].mData) + fixup;
		memcpy(mColumns[attribute].mData, a.defaultValue, t.size);
	}
}

inline FlexibleSOATable::FlexibleSOATable(const FlexibleSOATable& rhs) :
	mMeta(rhs.mMeta),
	mColumns(nullptr),
	mRowCapacity(0),
	mRowCount(0)
{
	assert(mMeta.GetAttributeCount() > 0);

	ReserveRows(rhs.mRowCapacity);
	for (AttributeID attribute = 0; attribute != mMeta.GetAttributeCount(); attribute++) {
		const MetaType& t = mMeta.GetType(mMeta.GetAttribute(attribute).type);
		memcpy(mColumns[attribute].mData, rhs.mColumns[attribute].mData, t.size * (rhs.mColumns[attribute].mUniform ? 1 : rhs.mRowCount));
		mColumns[attribute].mUniform = rhs.mColumns[attribute].mUniform;
	}
	mRowCapacity = rhs.mRowCapacity;
	mRowCount = rhs.mRowCount;
}

inline FlexibleSOATable::~FlexibleSOATable() {
	_aligned_free(mColumns[0].mData);
	free(mColumns);
}

inline const MetaTable& FlexibleSOATable::GetMetaTable() const {
	return mMeta;
}

inline size_t FlexibleSOATable::GetRowCapacity() const {
	return mRowCapacity;
}

inline size_t FlexibleSOATable::GetRowCount() const {
	return mRowCount;
}

inline void FlexibleSOATable::ReserveRows(size_t rowCapacity) {
	if (rowCapacity > mRowCapacity) {
		size_t newSize = 0;
		Column* columns = (Column*)malloc(mMeta.GetAttributeCount() * sizeof(Column));
		for (AttributeID attribute = 0; attribute != mMeta.GetAttributeCount(); attribute++) {
			const MetaType& t = mMeta.GetType(mMeta.GetAttribute(attribute).type);
			// newSize = (newSize + mMeta.GetAOSAlignment() - 1) & ~(mMeta.GetAOSAlignment() - 1);
			newSize = (newSize + t.alignment - 1) & ~(t.alignment - 1);
			columns[attribute].mData = (void*)newSize;	// Store offset first
			newSize += t.size * rowCapacity;
		}

		size_t fixup = (size_t)_aligned_malloc(newSize, mMeta.GetAOSAlignment());

		for (AttributeID attribute = 0; attribute != mMeta.GetAttributeCount(); attribute++) {
			const MetaType& t = mMeta.GetType(mMeta.GetAttribute(attribute).type);
			columns[attribute].mData = static_cast<char*>(columns[attribute].mData) + fixup;
			if (mColumns) {
				memcpy(columns[attribute].mData, mColumns[attribute].mData, t.size * (IsColumnUniform(attribute) ? 1 : mRowCount));
				columns[attribute].mUniform = mColumns[attribute].mUniform;
			}
		}

		if (mColumns) {
			_aligned_free(mColumns[0].mData);
			free(mColumns);
		}

		mColumns = columns;
		mRowCapacity = rowCapacity;
	}
}

inline void FlexibleSOATable::ReserveRowsInternal(size_t rowCount) {
	if (mRowCapacity < rowCount) {
		size_t newCapacity = mRowCapacity == 0 ? 16 : mRowCapacity * 2;
		while (newCapacity < rowCount)
			newCapacity *= 2;
		ReserveRows(newCapacity);
	}
}

inline void FlexibleSOATable::AppendRows(size_t rowCount) {
	size_t newCount = mRowCount + rowCount;
	ReserveRowsInternal(rowCount);

	// Initialize new rows
	size_t start = mRowCount;
	mRowCount = newCount; // Need to do this before calling SetValue()
	for (AttributeID attribute = 0; attribute != mMeta.GetAttributeCount(); attribute++) {
		// Only initialize the rows to default values if the column is varying.
		if (!IsColumnUniform(attribute)) {
			const void* defaultValue = mMeta.GetAttribute(attribute).defaultValue;
			for (size_t rowIndex = start; rowIndex < newCount; rowIndex++)
				SetValueRaw(rowIndex, attribute, defaultValue);
		}
	}
}

inline void FlexibleSOATable::InsertRows(size_t rowIndex, size_t rowCount) {
	assert(rowIndex < mRowCount);

	// To be optimized: reallocate and move first/second half.
	const size_t newCount = mRowCount + rowCount;
	ReserveRowsInternal(newCount);

	// Move second half to the back.
	for (AttributeID attribute = 0; attribute != mMeta.GetAttributeCount(); attribute++) {
		// Only move the rows if the column is varying.
		if (!IsColumnUniform(attribute)) {
			const size_t size = mMeta.GetType(mMeta.GetAttribute(attribute).type).size;
			char* p = static_cast<char*>(mColumns[attribute].mData) + rowIndex * size;
			memmove(p + rowCount * size, p, (mRowCount - rowIndex) * size);
		}
	}

	// Initialize new rows
	const size_t start = rowIndex;
	const size_t end = rowIndex + rowCount;
	mRowCount = newCount; // Need to do this before calling SetValue()
	for (AttributeID attribute = 0; attribute != mMeta.GetAttributeCount(); attribute++) {
		// Only initialize the rows to default values if the column is varying.
		if (!IsColumnUniform(attribute)) {
			const void* defaultValue = mMeta.GetAttribute(attribute).defaultValue;
			for (size_t i = start; i != end; ++i)
				SetValueRaw(i, attribute, defaultValue);
		}
	}
}

inline void FlexibleSOATable::RemoveRows(size_t rowIndex, size_t rowCount) {
	assert(rowIndex + rowCount <= mRowCount);

	if (rowCount > 0) {
		// move second half
		for (AttributeID attribute = 0; attribute != mMeta.GetAttributeCount(); attribute++) {
			// Only move the rows if the column is varying.
			if (!IsColumnUniform(attribute)) {
				const size_t size = mMeta.GetType(mMeta.GetAttribute(attribute).type).size;
				char* p = static_cast<char*>(mColumns[attribute].mData) + rowIndex * size;
				memmove(p, p + rowCount * size, (mRowCount - rowIndex - rowCount) * size);
			}
		}

		mRowCount -= rowCount;
	}
}

inline const void* FlexibleSOATable::GetValueRaw(size_t rowIndex, AttributeID attribute) const {
	if (IsColumnUniform(attribute))
		return mColumns[attribute].mData;
	else
		return const_cast<FlexibleSOATable *>(this)->GetValueRaw(rowIndex, attribute);
}

inline void* FlexibleSOATable::GetValueRaw(size_t rowIndex, AttributeID attribute) {
	ConvertColumnToVarying(attribute);

	const MetaAttribute& a = mMeta.GetAttribute(attribute);
	const MetaType& t = mMeta.GetType(a.type);
	return static_cast<char*>(mColumns[attribute].mData) + rowIndex * t.size;
}

inline void FlexibleSOATable::SetValueRaw(size_t rowIndex, AttributeID attribute, const void* value) {
	ConvertColumnToVarying(attribute);

	const MetaAttribute& a = mMeta.GetAttribute(attribute);
	const MetaType& t = mMeta.GetType(a.type);
	memcpy(static_cast<char*>(mColumns[attribute].mData) + rowIndex * t.size, value, t.size);
}

template <typename T>
const T& FlexibleSOATable::GetValue(size_t rowIndex, AttributeID attribute) const {
	assert(sizeof(T) == mMeta.GetType(mMeta.GetAttribute(attribute).type).size);
	return static_cast<const T*>(GetValueRaw(rowIndex, attribute))[0];
}

template <typename T>
T& FlexibleSOATable::GetValue(size_t rowIndex, AttributeID attribute) {
	assert(sizeof(T) == mMeta.GetType(mMeta.GetAttribute(attribute).type).size);
	return static_cast<T*>(GetValueRaw(rowIndex, attribute))[0];
}

template <typename T>
void FlexibleSOATable::SetValue(size_t rowIndex, AttributeID attribute, const T& value) {
	assert(sizeof(T) == mMeta.GetType(mMeta.GetAttribute(attribute).type).size);
	SetValueRaw(rowIndex, attribute, &value);
}

inline size_t FlexibleSOATable::GetAttributeSize() const {
	size_t newSize = 0;
	for (AttributeID attribute = 0; attribute != mMeta.GetAttributeCount(); attribute++) {
		const MetaType& t = mMeta.GetType(mMeta.GetAttribute(attribute).type);
		// newSize = (newSize + mMeta.GetAOSAlignment() - 1) & ~(mMeta.GetAOSAlignment() - 1);
		newSize = (newSize + t.alignment - 1) & ~(t.alignment - 1);
		newSize += t.size * mRowCapacity;
	}
	return newSize;
}

inline size_t FlexibleSOATable::GetMemorySize() const {
	return GetAttributeSize() 
		+ sizeof(this) 
		+ mMeta.GetAttributeCount() * sizeof(Column);// mColumns
}

inline bool FlexibleSOATable::IsColumnUniform(AttributeID attribute) const {
	assert(attribute < mMeta.GetAttributeCount());
	return mColumns[attribute].mUniform;
}

inline const void* FlexibleSOATable::GetUniformValueRaw(AttributeID attribute) const {
	assert(IsColumnUniform(attribute));
	return mColumns[attribute].mData;
}

inline void FlexibleSOATable::SetUniformValueRaw(AttributeID attribute, const void* value) {
	const MetaType& t = mMeta.GetType(mMeta.GetAttribute(attribute).type);
	memcpy(mColumns[attribute].mData, value, t.size);
	mColumns[attribute].mUniform = true;
}

template<typename T>
const T& FlexibleSOATable::GetUniformValue(AttributeID attribute) const {
	assert(sizeof(T) == mMeta.GetType(mMeta.GetAttribute(attribute).type).size);
	return static_cast<const T*>(GetUniformValueRaw(attribute))[0];
}

template<typename T>
void FlexibleSOATable::SetUniformValue(AttributeID attribute, const T& value) {
	assert(sizeof(T) == mMeta.GetType(mMeta.GetAttribute(attribute).type).size);
	SetUniformValueRaw(attribute, &value);
}

inline void* FlexibleSOATable::GetVaryingValuesRaw(AttributeID attribute) {
	ConvertColumnToVarying(attribute);
	return mColumns[attribute].mData;
}

inline void FlexibleSOATable::SetVaryingValuesRaw(size_t rowIndex, size_t rowCount, AttributeID attribute, const void* values) {
	assert(rowIndex + rowCount < mRowCount);
	const MetaType& t = mMeta.GetType(mMeta.GetAttribute(attribute).type);
	memcpy(static_cast<char*>(mColumns[attribute].mData) + t.size * rowIndex, values, t.size * rowCount);
	mColumns[attribute].mUniform = false;
}

inline void FlexibleSOATable::ConvertColumnToVarying(AttributeID attribute) {
	if (IsColumnUniform(attribute)) {
		// Copy uniform value to varying values
		const MetaType& t = mMeta.GetType(mMeta.GetAttribute(attribute).type);
		char* src = static_cast<char*>(mColumns[attribute].mData);
		char* end = src + t.size * mRowCount;
		for (char* dst = src + t.size; dst != end; dst += t.size)
			memcpy(dst, src, t.size);

		mColumns[attribute].mUniform = false;
	}
}