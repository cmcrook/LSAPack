#pragma once

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <string>
#include <thread>
#include <sstream>
#include <iostream>

namespace FG{
		using FloatingPointMicroseconds = std::chrono::duration<double, std::micro>;

		struct ProfileResult
		{
			std::string Name;
			FloatingPointMicroseconds Start;
			std::chrono::microseconds ElapsedTime;
			std::thread::id ThreadID;
		};

		struct ProfileSession
		{
			std::string Name;
		};

		class Profiler
		{
		public:
			inline Profiler(const Profiler&) = delete;
			inline Profiler(Profiler&&) = delete;

			inline void BeginSession(const std::string& name, const std::string& filepath = "results.json")
			{
				if (m_CurrentSession)
				{
					// If there is already a current session, then close it before beginning new one.
					// Subsequent profiling output meant for the original session will end up in the
					// newly opened session instead.  That's better than having badly formatted
					// profiling output.
					InternalEndSession();
				}
				m_OutputStream.open(filepath);

				if (m_OutputStream.is_open())
				{
					m_CurrentSession = new ProfileSession({ name });
					WriteHeader();
				}
				else
				{
					std::cout << "Profiler could not open results file!" << std::endl;
				}
			}

			inline void EndSession()
			{
				InternalEndSession();
			}

			inline void WriteProfile(const ProfileResult& result)
			{
				std::stringstream json;

				json << std::setprecision(3) << std::fixed;
				json << ",{";
				json << "\"cat\":\"function\",";
				json << "\"dur\":" << (result.ElapsedTime.count()) << ',';
				json << "\"name\":\"" << result.Name << "\",";
				json << "\"ph\":\"X\",";
				json << "\"pid\":0,";
				json << "\"tid\":" << result.ThreadID << ",";
				json << "\"ts\":" << result.Start.count();
				json << "}";

				if (m_CurrentSession)
				{
					m_OutputStream << json.str();
					m_OutputStream.flush();
				}
			}

			static Profiler& Get()
			{
				static Profiler instance;
				return instance;
			}
		private:
			inline Profiler()
				: m_CurrentSession(nullptr)
			{
			}

			inline ~Profiler()
			{
				EndSession();
			}

			inline void WriteHeader()
			{
				m_OutputStream << "{\"otherData\": {},\"traceEvents\":[{}";
				m_OutputStream.flush();
			}

			inline void WriteFooter()
			{
				m_OutputStream << "]}";
				m_OutputStream.flush();
			}

			// Note: you must already own lock on m_Mutex before
			// calling InternalEndSession()
			inline void InternalEndSession()
			{
				if (m_CurrentSession)
				{
					WriteFooter();
					m_OutputStream.close();
					delete m_CurrentSession;
					m_CurrentSession = nullptr;
				}
			}
		private:
			ProfileSession* m_CurrentSession;
			std::ofstream m_OutputStream;
		};

		class ProfileTimer
		{
		public:
			inline ProfileTimer(const char* name)
				: m_Name(name), m_Stopped(false)
			{
				m_StartTimepoint = std::chrono::steady_clock::now();
			}

			inline ~ProfileTimer()
			{
				if (!m_Stopped)
					Stop();
			}

			inline void Stop()
			{
				auto endTimepoint = std::chrono::steady_clock::now();
				auto highResStart = FloatingPointMicroseconds{ m_StartTimepoint.time_since_epoch() };
				auto elapsedTime = std::chrono::time_point_cast<std::chrono::microseconds>(endTimepoint).time_since_epoch() - std::chrono::time_point_cast<std::chrono::microseconds>(m_StartTimepoint).time_since_epoch();

				Profiler::Get().WriteProfile({ m_Name, highResStart, elapsedTime, std::this_thread::get_id() });

				m_Stopped = true;
			}
		private:
			const char* m_Name;
			std::chrono::time_point<std::chrono::steady_clock> m_StartTimepoint;
			bool m_Stopped;
		};

		namespace ProfileUtils {

			template <size_t N>
			struct ChangeResult
			{
				char Data[N];
			};

			template <size_t N, size_t K>
			constexpr auto CleanupOutputString(const char(&expr)[N], const char(&remove)[K])
			{
				ChangeResult<N> result = {};

				size_t srcIndex = 0;
				size_t dstIndex = 0;
				while (srcIndex < N)
				{
					size_t matchIndex = 0;
					while (matchIndex < K - 1 && srcIndex + matchIndex < N - 1 && expr[srcIndex + matchIndex] == remove[matchIndex])
						matchIndex++;
					if (matchIndex == K - 1)
						srcIndex += matchIndex;
					result.Data[dstIndex++] = expr[srcIndex] == '"' ? '\'' : expr[srcIndex];
					srcIndex++;
				}
				return result;
			}
		}

#define FG_PROFILE 1
#if FG_PROFILE
	// Resolve which function signature macro will be used. Note that this only
	// is resolved when the (pre)compiler starts, so the syntax highlighting
	// could mark the wrong one in your editor!
#if defined(__GNUC__) || (defined(__MWERKS__) && (__MWERKS__ >= 0x3000)) || (defined(__ICC) && (__ICC >= 600)) || defined(__ghs__)
#define FG_FUNC_SIG __PRETTY_FUNCTION__
#elif defined(__DMC__) && (__DMC__ >= 0x810)
#define FG_FUNC_SIG __PRETTY_FUNCTION__
#elif (defined(__FUNCSIG__) || (_MSC_VER))
#define FG_FUNC_SIG __FUNCSIG__
#elif (defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 600)) || (defined(__IBMCPP__) && (__IBMCPP__ >= 500))
#define FG_FUNC_SIG __FUNCTION__
#elif defined(__BORLANDC__) && (__BORLANDC__ >= 0x550)
#define FG_FUNC_SIG __FUNC__
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901)
#define FG_FUNC_SIG __func__
#elif defined(__cplusplus) && (__cplusplus >= 201103)
#define FG_FUNC_SIG __func__
#else
#define FG_FUNC_SIG "HZ_FUNC_SIG unknown!"
#endif

#define FG_PROFILE_BEGIN_SESSION(name, filepath) ::FG::Profiler::Get().BeginSession(name, filepath)
#define FG_PROFILE_END_SESSION() ::FG::Profiler::Get().EndSession()
#define FG_PROFILE_SCOPE_LINE2(name, line) constexpr auto fixedName##line = ::FG::ProfileUtils::CleanupOutputString(name, "__cdecl ");\
											   ::FG::ProfileTimer timer##line(fixedName##line.Data)
#define FG_PROFILE_SCOPE_LINE(name, line) FG_PROFILE_SCOPE_LINE2(name, line)
#define FG_PROFILE_SCOPE(name) FG_PROFILE_SCOPE_LINE(name, __LINE__)
#define FG_PROFILE_FUNCTION() FG_PROFILE_SCOPE(HZ_FUNC_SIG)
#else
#define FG_PROFILE_BEGIN_SESSION(name, filepath)
#define FG_PROFILE_END_SESSION()
#define FG_PROFILE_SCOPE(name)
#define FG_PROFILE_FUNCTION()
#endif

}